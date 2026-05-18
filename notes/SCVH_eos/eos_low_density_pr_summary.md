# EOS low-density PR summary

This is the short PR-facing summary for the low-density EOS work. The longer
working log is in `notes/SCVH_eos/mesa_scvh_low_pressure_extension.md`.

## What changed

The packaged SCVH input tables are left at their restored upstream extent. The
low-pressure SCVH coverage is now supplied analytically in
`eos/eosDT_builder/src/scvh_core.f90`, below the processed table edge at
`logP_min = -0.6`.

For each fixed `logT` row and each species, the tail uses the edge state

```text
logP0 = -0.6
logrho0 = logrho(logP0, logT)
logu0 = logu(logP0, logT)
logs0 = logs(logP0, logT)
R = 10**(logP0 - logrho0 - logT)
```

and the ideal low-pressure limit

```text
logrho_ideal(logP) = logrho0 + logP - logP0
logu_ideal(logP) = logu0
S_ideal(logP) = S0 + R*ln(10)*(logP0 - logP)
logs_ideal(logP) = log10(S_ideal)
comp_ideal(logP) = comp0
```

The implementation uses this ideal value below the first lower-pressure grid
step and applies a local cubic correction between `logP0` and that point so the
tail is value-continuous and has the same `d/dlogP` as the restored processed
SCVH table at the edge. That keeps the generated eosDT interpolation support
smooth without editing the source SCVH data files.

The regenerated OPAL/SCVH eosDT support no longer falls back to HELM for finite
SCVH pressure-ionization structure. That fallback produced noisy `mu` features
near the OPAL/SCVH support boundary in regenerated tables. In
`eos/eosDT_builder/src/helm_opal_scvh_driver.f90`, OPAL/SCVH points are rejected
only for non-finite or unphysical support values, and the last-resort rectangular
support fill uses an ideal HELM-compatible value when HELM itself is outside its
coverage.

FreeEOS fallback rows were corrected in
`eos/eosFreeEOS_builder/src/free_eos_table.f90`. Rows generated from the MESA EOS
store `lnPgas`, `lnE`, and `lnS` in natural-log form, but the FreeEOS data files
store those columns as base-10 logarithms. The fixed fallback path writes

```text
log10(Pgas) = lnPgas/ln10
log10(E) = lnE/ln10
log10(S) = lnS/ln10
log10(rho) = log10Rho
```

This removes the bad low-density `lnE` partial stripe near the
`logT = 3.0..3.1` FreeEOS boundary.

## Data and versions

The PR includes regenerated EOS data archives:

- `eos/eosFreeEOS_data.tar.xz`
- `eos/eosDT_data.tar.xz`

The FreeEOS table version and the eosDT table version are both bumped to `52`.
This is intentional: regenerated tables must invalidate existing version-51
binary caches.

## Region changes

The eosDT region definitions in `eos/eosDT_builder/src/eos_regions_defs.dek` and
`eos/eosDT_builder/src/eos_regions_code.dek` keep the generated OPAL/SCVH support
out of the warm, very-low-density FreeEOS side. This prevents the SCVH
low-pressure floor from interacting with the FreeEOS boundary at
`logrho ~ -15`.

## Verification notes

The current diagnostic images are kept in `notes/SCVH_eos/`:

- `eos_plotter_mu_current_v52_freeeos_fix.png`
- `eos_plotter_lnE_dlnT_v52_freeeos_fix.png`
- `eos_plotter_frac_OPAL_SCVH_current.png`
- `eos_plotter_frac_HELM_current.png`

The FreeEOS unit fix was checked against fallback rows near `logrho = -15`:

```text
logT = 3.02  FreeEOS row  logE = 12.9332
logT = 3.00  fallback row logE = 11.9080
logT = 2.98  fallback row logE = 11.8322
```

The remaining high-temperature, very-low-density region is not an SCVH target
region. The updated eosDT region mask leaves that side to FreeEOS/ideal coverage
instead of forcing OPAL/SCVH or HELM support there.
