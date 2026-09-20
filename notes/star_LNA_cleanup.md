# Star LNA Cleanup

Status: rebased onto the hydro changes branch. Source integration, static
checks, and a fresh MESA install are complete. MESA model validation has not
been performed.

Branch: `EbF/star_lna`

Rebase baseline: `origin/EbF/hydro_changes_and_bug_fixes` at `0a14f63c`,
including face reconstruction, the harmonic TDC dissipation length, and the
hydro energy-conservation changes.

## Scope

- [x] Update branch and implementation status in the Star LNA notes.
- [x] Add an explicit problem object and equation registry; move map and matrix
  lifetime management into problem setup and cleanup; report equation counts.
- [x] Extend row audits for RSP2 terms and TDC face placement.
- [x] Clarify growth and work output conventions; use the operator EOS and RSP2
  pressure closures in the work output.
- [x] Add a stop after LNA control, document the MLT treatment, and give complete
  saved model workflows.
- [x] Add a LaTeX manuscript deriving the implemented equations and mapping each
  equation to its source routine.
- [x] Complete static source checks and review the final branch diff against
  `origin/main`.

## Excluded Work

MESA model runs remain tasks for the user. Quantitative validation against RSP/RSP2
work files and nonlinear kicked models is not part of this cleanup pass.

## Face Reconstruction Invariants

For `use_face_reconstruction = .true.`, the Star LNA TDC closure must use the
same reconstructed face EOS and opacity state as the nonlinear MLT/TDC path.
For `use_face_reconstruction = .false.`, it must retain the stored face state
path. In both cases the face quantities must remain automatic differentiation
objects so their derivatives enter the linear operator.

The zone convention is `k = 1` at the surface and `k = nz` at the inner
boundary. Face audits must state the zone and face index explicitly.

## Work Log

- Generic pulsation case, 2026-08-27: renamed the tracked
  `dev_TDC_Cepheid_Hertzsprung_progression` case to
  `dev_TDC_Cepheid_Pulsation`. The new case uses Star LNA rather than embedded
  GYRE, includes the 6 M helium-depletion model, exposes both hydrodynamic
  velocity flags in `&star_job`, and lists the standard and `dPrad/dm`
  temperature-gradient forms. `x_integer_ctrl(8)` now sets the requested
  timesteps per period, replacing the fixed divisor of 600 in the local
  extras source. Star LNA initial velocity setup now writes either the `u` or
  `v` hydro state and rebuilds the Riemann face state for `u_flag`.
- The generic pulsation inlist activates
  `harmonic_dissipation_length_beta = 1`, so nonlinear TDC and Star LNA use
  the same harmonic length,
  \(1/\Lambda=1/(\alpha_{\rm MLT}H_P)+1/(\beta r)\). Its final `u_flag`
  section carries the PPISN metric split/merge AMR scheme and conservative
  child-pressure reconstruction. The supplied case freezes the mesh, so the
  scheme remains inactive until both `use_split_merge_amr` and
  `okay_to_remesh` are enabled.
- The generic inlist sets `min_kap_for_dPrad_dm_eqn = 1d-4` with the
  temperature-gradient controls. This opacity floor applies when the dPrad/dm
  temperature-gradient equation is enabled.
- The generic inlist also enables `floor_momentum_outer_BC_at_Prad`. With
  `use_momentum_outer_BC = .true.`, the surface boundary applies
  \[
  P_{\rm bc}\leftarrow\max\left(P_{\rm bc},\frac{aT_{\rm bc}^4}{3}\right).
  \]
  This is the PPISN surface-pressure floor that prevents a negative implied
  gas pressure. It is independent of the dPrad/dm opacity floor.
- `convergence_ignore_equL_residuals` and
  `make_gradr_sticky_in_solver_iters` are listed as temperature-gradient
  options. The former is commented out. `x_logical_ctrl(25)` retains the
  startup override in `run_star_extras`: the equL residual is ignored for the
  first ten models and included afterward. With the control false, the extras
  source leaves the native convergence flag unchanged.
- Hydro rebase, 2026-08-27: rebased the three Star LNA commits onto
  `0a14f63c`. The rebased commits are `19752ff4`, `071ea1a5`, and `08b7ba67`.
  The safety branch `EbF/star_lna_pre_hydro_rebase_20260827` retains the old
  tip `c86d36d3`.
- Hydro compatibility, 2026-08-27: `constant_L` remains a hydro-only equation.
  `check_star_LNA_model` rejects it before matrix assembly. Hydro
  `eval_dwork`, hydro `eval_dlnPdm_qhse`, and TDC `compute_Chi_cell` remain
  private. Star LNA uses its own static pressure-gradient and eddy-viscosity
  forms.
- TDC length integration, 2026-08-27: Star LNA now obtains the TDC mixing
  length from `get_mlt_mixing_length`. The source and convective luminosity use
  the effective `Lambda/Hp`, while damping, radiative loss, and eddy viscosity
  use `Lambda`. This matches the nonlinear TDC closure when the harmonic
  dissipation length is active.
- GYRE length policy, 2026-08-27: schema 130 still exports the legacy
  `mixing_length_alpha` and `Hp_face` fields. The writer warns when harmonic
  dissipation is active because GYRE will reconstruct `alpha*Hp`, not MESA's
  harmonic `Lambda`.
- Hydro energy audit, 2026-08-27: the hydro branch adds total-energy work
  proportional to `v*Uq` together with `Eq`. On a static radial background,
  `v`, `Uq`, and the velocity gradient are first order. Both `v*Uq` and `Eq`
  are therefore second order, so the current zero first-order Star LNA heating
  source remains correct.
- Main merge and dev case review, 2026-08-15: merged `origin/main` at
  `ee1f16e2`. Upstream removed the alternate density form TDC eddy viscosity
  control and implementation. The merge accepts that removal and deletes the
  obsolete Star LNA documentation and setup message for the old control.
- The 6M and 9M cases from main remain unchanged. A separate untracked
  `dev_TDC_Cepheid_6M_LNA` baseline now lives under
  `Classical_Pulsations/star_lna/models/star/dev_cases_TDC_Pulsation/`. It uses
  the updated 6M remesh, runs Star LNA at model 200, writes 15 modes, applies a
  5 km/s fundamental mode kick, and continues the nonlinear TDC calculation.
  The GYRE kick is disabled, while periodic GYRE analysis remains available.
- Static review found a duplicate `read_extra_pgstar_inlist(1)` assignment in
  the upstream 9M pulse header. The second assignment replaces the first, so
  the `pgstar` block in `inlist_pulses` is not read. This is left unchanged in
  the Star LNA branch because the shared 9M case is outside this branch's
  scope.
- RSP2 cell scale height pass, 2026-08-13: replaced the average of face pressure
  scale heights with `Peos(k)/(rho(k)*grav_cell)`, where `grav_cell` is the
  average gravity at the bounding faces. The nonlinear residual, automatic
  differentiation Jacobian, pre-solver turbulent velocity estimate, eddy
  viscosity coefficient, and Star LNA operator use the same definition.
- Documentation punctuation pass, 2026-08-12: removed decorative dash
  separators and unnecessary compound hyphens from the Star LNA prose. Minus
  signs in equations and literal command names were not changed.
- Documentation pass, 2026-08-12: revised the Star LNA and schema 130 control
  descriptions to match nearby MESA defaults. Restrictions, units, equations,
  output filenames, and mode indices are stated directly. Runtime messages,
  output headers, and the mathematical manuscript use the same terms. The
  detailed source cleanup note retains the original row designs, equations,
  source comparisons, edit sequence, and implementation status. The
  implementation note corrects the velocity kick limit to three selected modes.
  All note files remain untracked in the MESA checkout.
- Rebasing onto `4e05d577` preserved the Star LNA patch. The only textual
  conflict was the development controls declaration block; both upstream and
  Star LNA controls were retained.
- The matrix summary reports row norms and background residuals. The row,
  RSP2-term, and TDC face audits use the same output prefix.
- `star_LNA_problem` now owns the variable map, equation map, and dense matrix.
  The matrix summary reports equation counts, and the row structure output names
  the equation occupying every matrix row.
- RSP2 matrix summary runs now write `<prefix>_rsp2_term_audit.data`. TDC runs
  write `<prefix>_tdc_face_audit.data`, including active reconstructed and stored
  face states for direct comparison.
- TDC thermodynamic and opacity inputs pass through
  `get_reconstructed_face_state_ad`. With `use_face_reconstruction = .true.`,
  the LNA uses the reconstructed face EOS and opacity state.
- `star_LNA_stop_after_run` supports load model, analyze, and stop workflows in
  the normal and multiple star drivers. The binary startup path rejects the control
  explicitly because it cannot consume the stop result. The default remains
  false.
- `star_LNA_convection_treatment = 'mlt_static'` names the existing static MLT
  temperature gradient closure and is rejected for TDC and RSP2 backgrounds.
- Mode selection and output now name the filtered quantity
  `logKE_per_cycle`. The legacy `*_eta` control names are retained for inlist
  compatibility.
- Pressure and RSP2 turbulent pressure perturbations are built by the same
  closure for operator and work output. TDC MLT turbulent pressure remains a
  distinct face centered momentum term.

## Validation Log

- A fresh MESA install completed with exit status 0 on 2026-08-27 after
  restoring the Git LFS objects skipped during the rebase. The install built
  the Star LNA modules and passed the `star`, `astero`, and `binary` checks. No
  MESA model was run.
- The fresh build found one missing `wrap_r_p1` import in
  `star_LNA_turbulence_closures.f90`. The import was added and the resumed
  install passed.
- The 2026-08-27 control audit found every `star_LNA_*` control and
  `gyre_write_tdc_lna_background` in the development defaults, `star_info`
  storage, namelist declaration, input assignment, and output assignment.
- The hydro rebase retains the RSP2 face-temperature correction, cell pressure
  scale height, `w^2` radiative damping, `RSP2_nz_div_IBOTOM` cutoff, and
  turbulent-pressure normalization corrections.
- `git diff --check`, the added-line 132-column check, and Fortitude passed for
  the rebased source on 2026-08-27.

- Merge commit `c86d36d3` has no unmerged index entries. `git diff --check`
  passed, and Fortitude passed all 13 branch-specific `.f90` files relative to
  `origin/main`. Include fragments were not passed as standalone Fortran
  translation units.
- The external `dev_TDC_Cepheid_6M_LNA/standard_he_dep.mod` checksum is
  `a84daed587eca4e68a8dc2bb2217eaf2cb9eb19c1259f299edcf39542e6bd4d0`,
  identical to the updated 6M baseline from main. The external copy matches the
  reviewed temporary source tree.
- The generic `dev_TDC_Cepheid_Pulsation` case also includes the 80 M
  `model_00000800.mod` used by the star-top comparison. Its SHA-256 checksum is
  `c86e0e3c7a03c4de55f81f3b829de38d91a49dc82b1bf6d0775c460da10ef079`.
  The 6 M model remains active, and the 80 M load filename, GS98 opacity
  choices, `Zbase = 0.02`, and 5 MK center removal are documented as the
  alternative configuration.
- A fresh install on 2026-08-27 used
  `MESA_DIR=/Users/owner/Documents/Software/dev/test/mesa` and completed with
  exit status 0. The `star`, `astero`, and `binary` checker outputs matched
  their references. No stellar model was run.
- The RSP2 cell scale height change passes `fortitude check`, `git diff
  --check`, and the 132 column source check. The updated manuscript builds
  without LaTeX warnings, and the revised RSP2 closure page was inspected.
- GYRE interface audit, 2026-08-12: compared the MESA schema 130 writer with
  `EbF/gyre_convection_TDC` at `c61f5977`. The local checkout, local remote
  reference, and remote branch head all identify this commit.
- MESA and GYRE both set the schema 130 point data count to 38. MESA columns
  20 through 38 are read by GYRE in the same order: `L_conv0`, `A0`, `D_h`,
  `Hp`, `alpha_mlt`, `Cp`, `chiT`, `chiRho`, `gradL`, `gradT`, `alpha_C`,
  `alpha_S`, `alpha_D`, `alpha_R`, `alpha_Pt`, `alpha_M`, the MLT correction
  flag, `mlt_Pturb_factor`, and `mlt_Pturb0`.
- The exported velocity and turbulent pressure use the conventions required by
  GYRE,
  \[
  A_0 = v_{\rm conv}/\sqrt{2/3}, \qquad
  P_{{\rm turb},0} = f_{\rm Pt}\rho v_{\rm conv}^2/3.
  \]
  The exported horizontal diffusivity follows
  \[
  D_{\rm conv,h} =
  \frac{L_{\rm conv}H_P}
       {4\pi r^2\rho T c_P(\nabla_T-\nabla_L)},
  \]
  in cgs units. Atmosphere and center rows retain the repeated controls and set
  local TDC background values to zero, as expected by the GYRE reader.
- `git range-diff` confirmed that the rebased Star LNA commit differs from its
  before rebase form only where upstream added development controls.
- `fortitude check` passed for the changed Star LNA and driver Fortran files.
- `git diff --check` and the 132-column Fortran line length check passed after
  the source edits.
- The edited source and note files contain no en or em dash characters. The
  rendered manuscript text contains none.
- `star_LNA_manuscript.tex` builds without LaTeX warnings. The 10-page PDF was
  checked with `pdfinfo`, extracted with `pdftotext`, and visually inspected.
- The final branch diff against `origin/main` contains Star LNA and the
  schema 130 GYRE background export. It contains no envelope builder files and
  no files under `notes/`.
- The untracked Star LNA notes and manuscript are copied to
  `Classical_Pulsations/star_lna/notes/`. The untracked GYRE Markdown notes are
  copied to its `gyre/` subdirectory.
- On 2026-08-12, `make star` completed with exit status 0 using GNU Fortran
  15.2.0. `star_LNA_turbulence_closures.f90`, `star_LNA_support.f90`, and
  `star_LNA.f90` all compiled, and `build/star/lib/libstar.a` was rebuilt. The
  compiler emitted warnings in existing MESA sources, but none in the three
  Star LNA modules.
- No MESA model run has been performed.
- A user run of `dev_TDC_Cepheid_Pulsation` on 2026-08-27 loaded the 6 M
  helium-depleted model, removed the center at 2 MK, and completed the
  automatic 190-zone TDC remesh. It then stopped before the first step because
  the renamed case did not register `extras_check_model`.
- The complete TDC pulsation hook set has been restored in
  `run_star_extras.f90`: check-model, after-evolve, history-column,
  profile-column, startup, start-step, finish-step, and photo hooks are all
  registered. The corresponding no-op TDC after-evolve and profile-column
  interfaces were restored without reintroducing GYRE initialization,
  finalization, controls, or velocity kicks.
- `fortitude check`, `git diff --check`, and the 132-column source check pass
  after restoring the hooks. The assistant did not rebuild or rerun the work
  directory.
- The 2026-08-27 dPrad/dm run did not accept model 1 before its output ended.
  Solver call 1 used the requested initial timestep of (10^7\) s and failed
  after an excessive temperature correction. Calls 2 through 8 retried with
  timesteps from (10^5\) s down to (1.5625\times10^3\) s. The printed
  maximum residual was `dv_dt` or `dlnE_dt`, not `equL`. Call 8 reached a
  maximum residual of (3.21\times10^{-6}\), below the (10^{-5}\) maximum
  tolerance, but its residual norm was (1.63\times10^{-8}\), slightly above
  the (10^{-8}\) tolerance. The trace ended during call 9.
- The startup problem is a coupled structural relaxation after core removal
  and the 190-zone remesh. The dPrad/dm row remains in the Newton matrix even
  while `convergence_ignore_equL_residuals` omits it from the convergence
  norms. The trace does not show `equL` as the residual that rejects the later
  attempts.
- `x_ctrl(18)` was raised from (2\times10^3\) s to (10^6\) s. After model
  `x_ctrl(13)`, the effective pulsation limit is the smaller of this fallback
  and `dynamic_timescale/x_integer_ctrl(8)` or
  `period/x_integer_ctrl(8)`. This change does not affect model 1.
- `restore_mesh_on_retry` and `num_steps_to_hold_mesh_after_retry` were moved
  from the `u_flag` block to the TDC remesh section. The current inlist
  disables standard mesh adjustment after the one-time TDC remesh. The retry
  controls remain available when standard mesh adjustment is enabled.
- The logical `remesh_for_TDC_pulsations_when_load` control was replaced by
  the integer `steps_before_remesh_for_TDC_pulsations`. A negative value
  disables the remesh, zero retains the pre-first-step behavior in
  `run_star_support.f90`, and a positive value defers the remesh in
  `prepare_for_new_step` until
  \[
  n_{\rm model}=n_{\rm init}+N_{\rm remesh}.
  \]
  Here `N_remesh` counts accepted steps after the loaded model. The deferred
  remesh runs before the retry mesh snapshot, so retries retain the new mesh.
  Standard mesh adjustment is skipped on the remesh step. Photos preserve
  `init_model_number`, so a restart before the target retains the pending
  remesh and a restart after the target does not repeat it.
- A full `./install` after the remesh timing change completed successfully on
  2026-08-27 with
  `MESA_DIR=/Users/owner/Documents/Software/dev/test/mesa` and MESA SDK GNU
  Fortran 15.2.0. The updated controls, `prepare_for_new_step`, and startup
  remesh path compiled. The `star`, `astero`, and `binary` checker outputs
  matched their references. No development-case model was run.
