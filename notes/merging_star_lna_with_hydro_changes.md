# Merging Star LNA With Hydro Changes

Status: complete.

Branch: `EbF/star_lna`

Target: `origin/EbF/hydro_changes_and_bug_fixes`

## Scope

Rebase the three Star LNA commits onto the hydro changes branch while keeping
the Star LNA implementation limited to static radial LNA models. Hydro-only
equations and development test cases remain inherited base behavior.

## Invariants

- `star_LNA` does not support `constant_L`; reject that option before matrix
  assembly instead of adding a constant-luminosity LNA row.
- Hydro's `eval_dwork`, `eval_dlnPdm_qhse`, and TDC `compute_Chi_cell` remain
  private implementation routines. Star LNA uses its own static linear forms.
- The alternate density-form TDC eddy-viscosity control and implementation
  remain removed.
- For TDC, the active mixing and dissipation length is

  \[
  \Lambda = \alpha_{\rm MLT} H_P
  \]

  when `harmonic_dissipation_length_beta <= 0`. For a positive harmonic
  parameter,

  \[
  \frac{1}{\Lambda}
  =
  \frac{1}{\alpha_{\rm MLT}H_{P,\rm hse}}
  +
  \frac{1}{\beta_{\rm harm}r}.
  \]

  Star LNA must use the same `Lambda` as nonlinear TDC. In particular, the
  eddy-viscosity coefficient contains `TDC_alpha_M*Lambda`.
- GYRE schema 130 continues to export the legacy pair
  `mixing_length_alpha`, `Hp_face`. If harmonic dissipation is active, emit one
  warning that GYRE reconstructs `Lambda = mixing_length_alpha*Hp_face`, which
  does not match the active MESA harmonic length.
- Preserve the existing RSP2 cell pressure-scale-height, radiative-damping,
  cutoff, and turbulent-pressure corrections.
- Keep `notes/` and generated output untracked.
- For `use_dPrad_dm_form_of_T_gradient_eqn`, use the same interior algebraic
  relation as `hydro_temperature`:

  \[
  \Delta P_{\rm rad,expected}
  =
  -\frac{\Delta m\,\kappa_f L_{\rm rad}}
  {c A_f^2\lambda_f},
  \qquad
  \Delta P_{\rm rad,actual}
  =
  \frac{a}{3}\left(T_{k-1}^4-T_k^4\right),
  \]

  where `lambda_f = 1` without flux limiting. Preserve the opacity floor,
  face reconstruction, MLT radiative-luminosity split, and flux limiter.
- For `u_flag`, the dynamic velocity variable is cell-centered `u_k`. Use
  the Riemann momentum equation

  \[
  \frac{d u_k}{dt}
  =
  \frac{F_{P,\mathrm{in}}-F_{P,\mathrm{out}}
  +F_{\mathrm{geometry}}+F_{\mathrm{gravity}}+F_{U_q}}
  {\Delta m_k},
  \]

  and the kinematic relation `d lnR_k/dt = u_face,k/r_k`. Reuse the normal
  hydro Riemann face reconstruction and momentum right-hand side. Apply the
  optional initial velocity kick to the active `u` or `v` hydro variable.

## Checklist

- [x] Fetch the current hydro branch.
- [x] Record the starting commits and create a local safety branch.
- [x] Rebase the Star LNA commits onto the hydro branch.
- [x] Keep the hydro removal of `TDC_use_density_form_for_eddy_viscosity`.
- [x] Remove obsolete density-form Star LNA comments and runtime references.
- [x] Remove unnecessary public exports for `eval_dwork`,
  `eval_dlnPdm_qhse`, and TDC `compute_Chi_cell`.
- [x] Preserve hydro's `do1_constant_L_eqn` without adding Star LNA support.
- [x] Reject `constant_L` in `check_star_LNA_model`.
- [x] Adapt `set_TDC_LNA` to the hydro branch's AD mixing-length interface.
- [x] Use the hydro mixing length in the Star LNA TDC eddy-viscosity closure.
- [x] Keep GYRE schema 130 on the legacy `alpha*Hp` convention and add one
  warning for harmonic dissipation.
- [x] Verify all Star LNA controls are declared, read, stored, and written.
- [x] Verify the RSP2 fixes remain present after the rebase.
- [x] Review the hydro energy-conservation changes against the static LNA
  linearization.
- [x] Update the Star LNA notes and manuscript equations.
- [x] Copy the updated notes to the external Star LNA notes directory.
- [x] Run `git diff --check`, line-length checks, and Fortitude.
- [x] Perform a fresh MESA install.
- [x] Review the final history, tracked diff, and untracked files.
- [x] Add the `dPrad/dm` STAR LNA temperature-gradient row.
- [x] Add the cell-centered `u_flag` variable and momentum path.
- [x] Audit eigenfunction and work diagnostics for the active velocity grid.
- [x] Update controls documentation and Star LNA notes.
- [x] Copy the updated notes to the external Star LNA notes directory.
- [x] Run non-compiling style, line-length, and diff checks.
- [x] Rename `use_P_d_1_div_rho_form_of_work_when_time_centering_velocity`
  to `use_P_d_1_div_rho_form_of_work`.
- [x] Apply the simple work form independently of velocity time centering.
- [x] Omit TDC `v*Uq` from the ordinary total-energy source until the hydro
  matrix supports its second-neighbor derivatives.
- [x] Use the selected `dPrad/dm` transport relation inside the perturbed TDC
  luminosity closure.
- [x] Report the row location of the first otherwise selectable root rejected
  by the full-pencil residual limit.
- [x] Copy the TDC energy-work note to the external Star LNA notes directory.
- [x] Run non-compiling checks for the TDC energy-work change.

## Work Log

- 2026-08-27: created this checklist before changing branch history.
- 2026-08-27: fetched `origin`. The starting Star LNA tip is `c86d36d3`,
  the hydro target is `0a14f63c`, and their merge base is `ee1f16e2`.
  Created local safety branch `EbF/star_lna_pre_hydro_rebase_20260827` at
  `c86d36d3`.
- 2026-08-27: rebased the three Star LNA commits onto `0a14f63c`. The new
  commits are `19752ff4`, `071ea1a5`, and `08b7ba67`.
- 2026-08-27: kept `do1_constant_L_eqn` as hydro behavior and added an early
  `constant_L` rejection to `check_star_LNA_model`.
- 2026-08-27: removed stale density-form references and the unnecessary
  public exports for hydro `eval_dwork`, hydro `eval_dlnPdm_qhse`, and TDC
  `compute_Chi_cell`. The hydro routines remain available to their internal
  callers. Star LNA retains its separate static HSE-gradient routine.
- 2026-08-27: updated `set_TDC_LNA` to accept the effective AD mixing-length
  alpha and `Lambda`. Source and convective-flux terms use the effective
  alpha, while turbulent and radiative damping use `1/Lambda`, as in the
  nonlinear TDC path.
- 2026-08-27: changed the Star LNA TDC eddy-viscosity coefficient to use
  `TDC_alpha_M*Lambda_cell`, where `Lambda_cell` is the same face-average used
  by nonlinear TDC.
- 2026-08-27: kept GYRE schema 130 on its existing `alpha_MLT` and `Hp_face`
  fields. The writer now warns when harmonic dissipation is active because
  the exported GYRE closure then differs from the active MESA closure.
- 2026-08-27: the first fresh install attempt found that the rebased RSP2 cell
  gravity expression used `wrap_r_p1` without importing it from
  `auto_diff_support`. Added the missing import before resuming the install.
- 2026-08-27: `git diff --check` passed. No added Fortran line exceeds 132
  columns. Fortitude passed all 13 relevant Fortran files.
- 2026-08-27: restored the skipped Git LFS objects required by the installer.
  The fresh MESA install then completed with exit status 0, including the
  `star`, `astero`, and `binary` checks. No MESA model was run.
- 2026-08-27: rebuilt the 10-page Star LNA manuscript without LaTeX warnings.
  Rendered and inspected all pages; no clipping, overlap, or broken equations
  were found.
- 2026-08-27: copied the updated notes, manuscript source, and manuscript PDF
  to `Classical_Pulsations/star_lna/notes/` and verified matching SHA-256
  checksums. The repository notes and generated output remain untracked.
- 2026-08-27: final history review confirms three Star LNA commits directly on
  `0a14f63c`. The branch is three commits ahead of the hydro target. The local
  tracked compatibility diff is limited to seven files.
- 2026-08-27: reopened the checklist to support the alternate `dPrad/dm`
  temperature-gradient equation and cell-centered `u_flag` hydrodynamics.
  The `u_flag` operator will reuse the Riemann face state and cell momentum
  right-hand side.
- 2026-08-27: added the `dPrad/dm` algebraic row with the same opacity floor,
  reconstructed face state, MLT radiative-luminosity split, and flux limiter
  as `hydro_temperature:do1_alt_dlnT_dm_eqn`.
- 2026-08-27: completed the `u_flag` variable map, radius row, Riemann cell
  momentum row, pressure-work stencil, AD mapping, eigenfunction output, work
  normalization, and background residual diagnostics. The Riemann pressure and
  geometry helpers accept an optional time-centering switch so STAR LNA uses
  the continuous equations while normal hydro behavior is unchanged.
- 2026-08-27: extended initial velocity setup to `u_flag`. The kick writes the
  cell velocity and its current and start solver state, then rebuilds the
  Riemann face velocity and copies it into `u_face_start`.
- 2026-08-27: added a static STAR LNA cell eddy-viscosity source for TDC
  `u_flag` models. RSP2 eddy viscosity with `u_flag` remains rejected because
  nonlinear RSP2 still inserts this term into the Riemann face velocity rather
  than the cell momentum source.
- 2026-08-27: updated the controls documentation, implementation plan,
  readiness checklist, compression history, and mathematical manuscript for
  the two supported paths. Rebuilt and inspected the 11-page manuscript PDF;
  LaTeX reported no warnings or bad boxes.
- 2026-08-27: Fortitude passed the five changed Fortran implementation files.
  `git diff --check` passed, and no changed source or controls line exceeds 132
  columns. These additions have not been compiled or exercised in a MESA run.
- 2026-08-27: copied the revised checklist, plans, manuscript source, and
  manuscript PDF to `Classical_Pulsations/star_lna/notes/`. Matching SHA-256
  checksums verified every copy.
- 2026-08-27: completed a MESA-style source review of the `u_flag` and
  `dPrad/dm` additions. Tightened the active turbulent-pressure guard, made
  STAR LNA's Riemann options explicit at the call sites, restored concise
  force and gravity comments, removed the author-tagged RSP2 TODO, and kept
  only helpers that share equation logic or remove repeated velocity-grid
  selection. Fortitude, line-length, and whitespace checks still pass.
- 2026-08-28: renamed the simple `P d(1/rho)` work control and made it
  independent of velocity time centering. The TDC ordinary total-energy path
  no longer adds `v*Uq` because its half-cell work requires second-neighbor
  derivatives which are not represented by `auto_diff_real_star_order1`.
  RSP2 behavior is unchanged. The derivation and remaining face-`Eq` work are
  recorded in `tdc_eddy_viscous_energy_work.md`.
- 2026-08-28: `git diff --check`, the controls/defaults linter, line-length
  checks, and Fortitude on the changed `.f90` files pass. A fresh `./install`
  completed successfully, including the `star`, `astero`, and `binary` checks.
  No MESA model was run.
- 2026-08-29: the 350-zone Cepheid TDC audit showed that perturbed TDC still
  used the QHSE spatial gradient while the model used `dPrad/dm`. The inner
  luminosity closure differed from the model luminosity by 43 percent.
  `actual_gradT_for_star_LNA` now obtains `Lrad` from the active `dPrad/dm`
  equation and passes `gradT = Lrad/L0` to `set_TDC_LNA`. RSP2 already uses
  the shared `dPrad/dm` residual with its reconstructed radiative luminosity.
  Added a residual-row diagnostic for the first rejected mode candidate.
- 2026-08-29: the non-`dPrad/dm` TDC path mixed interpolated pressure in the
  QHSE spatial gradient with reconstructed pressure in `L0`. The pressure
  mismatch at zone 298 was `2.476d-3`, consistent with the `2.392d-3`
  luminosity-row residual. `star_LNA_eval_dlnPdm_qhse` now uses reconstructed
  face pressure when face reconstruction is active, matching
  `hydro_temperature:eval_dlnPdm_qhse`.
- 2026-08-29: the corrected QHSE run reduced the background luminosity
  residual to `9.015d-9`, but retained 228 zero-`w` constraints in the dynamic
  generalized pencil. `partition_star_LNA_indices` now eliminates rows without
  time derivatives. The 350-zone model should therefore use 1172 dynamic
  variables and 928 algebraic variables while retaining all 122 active TDC
  `w` equations.
- 2026-08-29: a fresh MESA install completed successfully after the algebraic
  partition correction.
- 2026-08-29: the corrected model-200 eigensolve confirmed 1172 dynamic and
  928 algebraic variables for 122 active TDC `w` zones. It recovered the
  unstable 11.8313-day fundamental with `logKE_per_cycle = 8.209d-2`, but the
  reconstructed full-pencil residual remained `2.796d-2`. A later calculation
  with 118 active TDC zones reduced that residual to `4.395d-3` without
  changing the physical mode identification. Reduced-pencil and full-pencil
  residuals must be compared before changing the solver or acceptance limit.
- 2026-08-29: added full-pencil Newton refinement for candidate acoustic
  eigenpairs. Each correction uses one complex band factorization and two
  solves of `A - sigma*B`, with residual-decreasing backtracking. The refined
  root is retained only when it lowers the residual against the original
  unscaled equations. The Cepheid development inlist now requests `1d-8`.
