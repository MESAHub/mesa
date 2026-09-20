# star_LNA Validation Status

Base status: synchronized with `origin/EbF/star_lna` at `622075fbf` on 2026-09-18,
including `origin/main` at `fd396fd73`. The earlier source audit was performed
on 2026-09-16 at `f6d606939`; the surface-boundary correction is now included.
The run record below includes a successful 2026-08-29 install after the TDC
transport and algebraic-partition corrections. The corrected model-200 runs
recovered the 11.83-day fundamental but still exceeded the full-pencil
residual limit. Newton refinement and the later August 31/September 1 changes
are implemented; this record does not establish their model validation.
No MESA build or model run was performed for this documentation revision.

The manuscript presentation was revised on 2026-09-18 at the same source
commit. Its main equation sequence now pairs model equations with their
linearized rows: blue for dynamic rows and rust for algebraic rows. Closure
coefficients and diagnostics follow in appendices. This layout change does
not change the source or model-validation checks below.

A follow-up on 2026-09-18 restored the explicit `dedt`/`eps_grav` comparison,
local and flux pressure-work forms, discrete mechanical/turbulent inertia,
and the `u_flag`/`v_flag` distinctions. At that revision, LNA rejected the
`eps_grav` form. The 29-page PDF was rebuilt without LaTeX warnings, all
pages were visually checked, and the notes/output copies were synchronized.
After upstream synchronization on the same date, the boundary equations and
source status were updated and the 29-page PDF was rebuilt and checked again.
These are documentation checks, not new model validation.

## Independent-Y source revision, 2026-09-19

The current working tree replaces the RSP2 Hp unknown with independent signed
face Y and a forward flux row, sharing the public TDC scale/mixing lengths.
RSP2 `eps_grav` and cell-grid viscosity are now implemented. The manuscript
retains the detailed energy/work distinctions and describes these uncommitted
changes. TDC viscosity remains separate; the shared cell mixing length
was subsequently updated in both nonlinear hydro and LNA.
Source, syntax, and Python algebra checks pass. On 2026-09-19 the user
authorized `./install`: compilation and installer package checks passed with
`MESA_DIR=/Users/owner/Documents/Software/dev/test/mesa` and SDK 26.6.1.
Actual RSP2 AD partial checks, restart/remesh models, and stellar convergence
tests remain pending user execution. The current progress and user-run validation matrix are in
[rsp2_independent_Y_implementation.md](rsp2_independent_Y_implementation.md).
The equation-focused manuscript includes the signed RSP2 PII product rule,
shared cell-local mixing length, explicit dynamical-gradient response and the
nonlinear/continuous surface-boundary distinction. The current rebuild and
visual review are tracked in `rsp2_independent_Y_implementation.md`.
Remesh controls, case defaults, build history, and saved-model workflow
instructions remain in the development notes rather than the manuscript.
The notes/output copies match. The historical run checkboxes below retain
their original scope; no new MESA model or eigenproblem was run for this refresh.

## Base source checks, updated 2026-09-18

- [x] Match the manuscript to local commit `622075fbf` after upstream synchronization.
- [x] Document `star_LNA_T_inner` and the fixed perturbations below its cut.
- [x] Document full-pencil Newton refinement before mode acceptance.
- [x] Match RSP2 row selection, the signed PII response, and cubic damping.
- [x] Document `TDC_use_dynamical_gradL` and its Ledoux restriction.
- [x] Distinguish cell and face stresses for the two velocity grids.
- [x] Record the static MLT/RSP2 inner `dPrad/dm` spacing correction and the
  separate TDC helper, which still uses `s%dm_bar(k)`.
- [x] Correct schema-130 column 22 to diffusivity, with units `cm^2/s`.
- [ ] Validate the RSP2 gradient and signed PII derivatives on saved models.
- [ ] Validate `TDC_use_dynamical_gradL`, including its composition term.
- [ ] Compare eddy-viscosity acceleration and work for `v_flag` and `u_flag`.
- [ ] Check the optional inner TDC stress at `R_center > 0`.
- [ ] Check both transport helpers at the physical inner boundary and at a
  `star_LNA_T_inner` cut through a full model.
- [x] Integrate surface-boundary commit `39fd1207a` and document its equations.
- [ ] Validate the integrated surface-boundary correction on saved models.
- [ ] Derive and validate a discrete work balance against the eigenvalue growth;
  the arithmetic sum of the six diagnostic columns is not an established identity.

The remaining run checks below retain their recorded status. A source match
or a rebuilt PDF does not mark a model test complete.

## Saved Model Execution

`star_LNA` runs through the normal `star` executable after loading a MESA/star
saved model. No separate executable is required.

This path initializes the normal MESA/star driver, reads the inlist, loads the
model, and calls star LNA. The old RSP `do_LINA` path is not a general reader
for MESA/star `.mod` files; it runs after RSP setup has constructed an RSP model
in memory.

Use this `star_job` block for each saved model test:

```fortran
&star_job
   load_saved_model = .true.
   load_model_filename = 'model.mod'
/
```

With `star_LNA_model_number < 0`, the analysis runs after `extras_startup`.
With a nonnegative value, it runs after `extras_start_step` when
`s% model_number` matches and before that step is evolved.

`star_LNA_stop_after_run = .true.` skips evolution after a successful analysis
in the normal and multiple star drivers. Startup profile, initial model, and echo
actions still run. The binary driver rejects this startup stop option.

### MLT

Use an MLT model without TDC or RSP2:

```fortran
&controls
   star_LNA_flag = .true.
   star_LNA_model_number = -1
   star_LNA_stop_after_run = .true.
   star_LNA_set_initial_velocity = .false.
   star_LNA_convection_treatment = 'mlt_static'
/
```

### TDC

Retain the TDC controls used to create the saved model:

```fortran
&controls
   star_LNA_flag = .true.
   star_LNA_model_number = -1
   star_LNA_stop_after_run = .true.
   star_LNA_set_initial_velocity = .false.
   star_LNA_include_tdc = .true.
   star_LNA_convection_treatment = 'perturbed'
/
```

### RSP2

Retain the RSP2 controls used to create the saved model:

```fortran
&controls
   star_LNA_flag = .true.
   star_LNA_model_number = -1
   star_LNA_stop_after_run = .true.
   star_LNA_set_initial_velocity = .false.
   star_LNA_include_rsp2 = .true.
   star_LNA_convection_treatment = 'perturbed'
/
```

TDC and RSP2 tests require saved models whose convective state was produced by
the same convection controls. Enabling TDC or RSP2 only for the analysis does
not initialize that state.

Set `star_LNA_set_initial_velocity = .true.` only for a velocity initialization
test. This option writes the active `u` or `v` variable and its current and
start solver state. For `u_flag`, it also rebuilds the Riemann face velocity
and initializes `u_face_start`. The LNA solve itself does not require a hydro
velocity variable; velocity initialization requires either `u_flag` or
`v_flag`.

## Verified Source State

The build and lint results in this section are historical checks at the
recorded dates; they were not repeated against `f6d606939` in this audit.

- [x] `fortitude check` passes for the Star LNA and driver files.
- [x] `git diff --check` passes.
- [x] Changed Fortran and controls lines do not exceed 132 columns.
- [x] `make star` completed with exit status 0 on 2026-08-12 using GNU Fortran
  15.2.0.
- [x] A fresh MESA install completed with exit status 0 on 2026-08-27 after the
  hydro rebase.
- [x] A fresh MESA install completed with exit status 0 on 2026-08-28 after
  adding `star_LNA_max_eigenvector_residual`.
- [x] The problem object owns the variable map, equation map, and dense matrix.
- [x] Row structure output records the equation and dominant A/B entries for
  every matrix row.
- [x] RSP2 output includes `<prefix>_rsp2_term_audit.data`.
- [x] Active TDC output includes `<prefix>_tdc_face_audit.data`.
- [x] TDC uses the reconstructed face EOS and opacity state when
  `use_face_reconstruction = .true.`.
- [x] TDC uses the hydro mixing length for normal and harmonic dissipation.
- [x] Perturbed TDC obtains its spatial gradient from `dPrad/dm` when that
  temperature-gradient equation is active.
- [x] The ordinary QHSE spatial gradient uses reconstructed face pressure when
  face reconstruction is active.
- [x] RSP2 uses the shared `dPrad/dm` residual with `s% Lr_ad(k)`.
- [x] Residual rejection reports the first otherwise selectable root and its
  maximum-residual row location.
- [x] Every Star LNA and schema 130 control has a default, storage field,
  namelist entry, input assignment, and output assignment.
- [x] `constant_L` is rejected before matrix assembly.
- [x] `frozen_flux` omits the internal TDC `w` variable and uses
  `delta(Lrad + Lconv0 - L) = 0` in interior luminosity rows.
- [x] `mlt_static` selects the static temperature gradient row and rejects TDC
  and RSP2 backgrounds.
- [ ] MLT, TDC, and RSP2 saved model tests.
- [ ] Direct comparison with RSP or RSP2 periods, growth, work, and velocity
  initialization.

## First TDC Kick Test, 2026-08-28

The first continued TDC test does not validate the selected modes or the
velocity initialization. At model 200, the third selected root has a period of
1498.86 days and `logKE_per_cycle = 7.72888`. Its real displacement slice has
90 sign changes across the 372-zone envelope, including three sign changes
between zones 20 and 23 near `logT = 3.56`. The kick therefore changes from
about `+5.0` to `-5.3` km/s across the first of these adjacent zones. This is
not a clean second-overtone displacement profile.

The nonlinear response is consistent with that defect. The accepted timestep
falls from `2.18176d5` s at model 200 to `1.06525d2` s at model 201 after 11
retries. The imposed surface velocity is `0.903*csound`. The run then loses the
surface velocity within a few days, much less than one percent of the selected
period, so the response is a short surface transient rather than the selected
mode.

The written velocity eigenfunction also fails the radius-row identity

```text
delta_v/r = sigma*delta_lnR
```

on the static background. For the selected third root, the median amplitude
ratio between the two sides is about 103. All 15 selected roots have many
displacement sign changes and fail the same check. Period and growth output
from this solve must not be used until the eigensolution satisfies the assembled
rows.

The kick builder formerly replaced the complex velocity eigenfunction with a
signed displacement-amplitude profile. It now phase-aligns each selected mode
at the surface and combines the normalized velocity eigenfunctions directly.

Required solver checks:

- [x] Evaluate `A*x - sigma*B*x` for every candidate root before selection.
- [x] Scale the reduced pencil after algebraic elimination and undo its column
  scaling before reconstructing the algebraic variables.
- [x] Confirm after rebuilding that the radius-row identity is satisfied to
  the accepted componentwise tolerance.
- [ ] Identify radial order from a validated displacement eigenfunction instead
  of treating frequency-list position as radial order.
- [ ] Repeat the kick test at a linear amplitude after the eigenvector residual
  check passes.

The reproducible analysis is in
`star/dev_cases_TDC_Pulsation/dev_TDC_Cepheid_Pulsation/plotting`.

### Dense Solver Diagnosis, 2026-08-28

An isolated run used
`MESA_DIR=/Users/owner/Documents/Software/dev/test/mesa`. With the inlist value
`star_LNA_min_mode_frequency_uHz = 0.1`, the selector cannot include a period
longer than 115.74 days. The development inlist now uses

```text
star_LNA_min_mode_frequency_uHz = 2.314814814814815d-3
```

which permits periods through 5000 days. A 1000-day mode has a frequency of
about `0.01157` microHz and is inside this window.

Lowering the frequency floor did not repair the original solve. Its first
selected root had a period of 297.79 days and a median radius-row amplitude
mismatch of about 425. Its displacement also had 72 nodes above one percent of
the surface amplitude.

The full matrix has raw column norms spanning about 74 orders of magnitude.
The dense path scales the full pencil before forming the Schur complement

```text
Ared = Add - Ada*inv(Aaa)*Aad
Bred = Bdd - Bda*inv(Aaa)*Aad,
x_alg = -inv(Aaa)*Aad*x_dyn.
```

Algebraic elimination changes the row and column norms. The reduced pencil is
therefore equilibrated by alternating row and column max-norm scaling. Before
equilibration, `Bred` is divided by the model dynamical timescale, so `DGGEV`
solves for

```text
lambda = sigma*t_dyn.
```

If `C` is the accumulated reduced column scaling and `y` is the returned
eigenvector, reconstruction uses

```text
x_dyn = C*y
x_alg = -inv(Aaa)*Aad*x_dyn
```

and then reverses the full-pencil column scaling. Mode selection evaluates the
maximum componentwise backward residual against the original full pencil,

```text
max_i |(A*x - sigma*B*x)_i|
      / (sum_j |A_ij*x_j| + |sigma|*sum_j |B_ij*x_j|),
```

and rejects values above `star_LNA_max_eigenvector_residual`. This control is
an upper limit, not a floor. It defaults to `1d-3` and cannot disable the
validity check. Candidate modes with an initial componentwise residual at most
`1d-1`, or a balanced normwise residual at most `1d-6`, are first refined
against the scaled full pencil. This permits repair of tiny moment components
without weakening the final componentwise check. A final limit of `1d-1` is suitable
only for diagnostic output because it accepts a row defect as large as ten
percent. The selected-mode table records the final residual.
The implementation is in `star/private/star_LNA_support.f90`, in
`solve_dense_star_LNA`, `scale_star_LNA_matrix_for_solver`,
`scale_star_LNA_pencil`, `refine_star_LNA_eigenpair`, and
`star_LNA_eigenvector_residual_from_vector`.

For an approximate eigenpair, define

```math
r=(A-\sigma B)x.
```

One Newton correction satisfies

```math
(A-\sigma B)\,\delta x-Bx\,\delta\sigma=-r,
\qquad c^\dagger\delta x=0.
```

The code fixes the largest component of `x` and uses one complex band
factorization with two right-hand sides,

```math
y=(A-\sigma B)^{-1}(-r),\qquad
z=(A-\sigma B)^{-1}Bx,
```

```math
\delta\sigma=-\frac{y_j}{z_j},\qquad
\delta x=y+z\,\delta\sigma.
```

`ZGBTRF` and `ZGBTRS` operate on the original block-banded full pencil in its
scaled coordinates. Up to four accepted corrections are allowed. A short
backtracking search requires each correction to reduce the scaled residual,
and the refined eigenpair replaces the dense result only if it also reduces
the residual measured against the original unscaled equations.

The installed model-200 test selects a zero-node fundamental at 1690.53 days
with `logKE_per_cycle = 2.18005`, followed by a one-node mode at 503.91 days.
The fundamental has a maximum componentwise residual of `6.28d-5`; the median
radius-row residual is `1.65d-14`. The seven selected modes have maximum
residuals between `6.28d-5` and `6.29d-4`. A `1d-6` limit rejects these
validated modes and is too strict for the present unscaled full-pencil metric.
Roots above the default `1d-3` limit are excluded.

The corrected fundamental kick uses the complex velocity eigenfunction. A
requested 5 km/s surface velocity peaks at 5.23 km/s in the interior. The
surface value is `0.906*csound`, so this remains a near-sonic nonlinear kick.
A linear startup test should use a substantially smaller velocity before the
subsequent timestep behavior is used to judge the mode shape.

## Work Output

Each `<prefix>_work_<mode>.data` file contains

```text
pressure_work
turb_pressure_work
eddy_visc_work
rad_lum_work
conv_lum_work
turb_lum_work
total_work
```

`pressure_work` and `turb_pressure_work` use the same EOS and RSP2 pressure
closures as the operator. The file divides all six components by the modal
kinetic energy. `total_work` is their sum.

Required comparisons:

- [ ] Compare `pressure_work`, `turb_pressure_work`, and `eddy_visc_work` with
  RSP `LINA_work*.data`.
- [ ] Compare the modal kinetic energy normalization with RSP.
- [ ] Determine the sign convention for `rad_lum_work`, `conv_lum_work`, and
  `turb_lum_work` from matched RSP/RSP2 output.
- [ ] Add a separate face centered TDC MLT turbulent pressure work term if the
  matched output requires it.

Until these comparisons are complete, the sum of the six columns is not an independent
estimate of eigenvalue growth.

## Growth Convention

Selected mode output uses

```text
logKE_per_cycle = 4*pi*sigma_real/sigma_imag
KE_fractional_growth = exp(logKE_per_cycle) - 1
GREKM = 2*tanh(logKE_per_cycle/2)
amplitude_fractional_growth = exp(logKE_per_cycle/2) - 1
```

Required comparisons:

- [ ] Compare nonlinear radius amplitude growth with
  `amplitude_fractional_growth_per_period`.
- [ ] Compare symmetric nonlinear kinetic energy growth with
  `grekm_growth_per_period`.
- [ ] Exclude transient cycles after the velocity initialization.

## RSP2

The local implementation uses the selected temperature-gradient row and the
hydro AD face state. RSP2 uses the same signed PII in nonlinear hydro
and LNA; TDC uses its unsaturated LNA relation.

`<prefix>_rsp2_term_audit.data` writes `PII`, `Lc`, `Lt`, source, dissipation,
radiative damping, turbulent pressure work, turbulent luminosity divergence,
RHS, and inertia for each zone.

Required checks:

- [ ] Compare `Lt_ad` in the total energy and turbulent energy rows with
  `hydro_rsp2.f90`.
- [x] Form RSP2 cell `Hp` from cell pressure and density with gravity averaged
  from the bounding faces in hydro, pre-solver initialization, and star LNA.
- [ ] Compare radiative damping values with matched RSP and RSP2 models.
- [ ] Check the forced nonturbulent cutoffs for `w`, `Lc`, `Lt`, and source
  terms.
- [ ] Inspect selected mode `w` and `Hp` amplitudes and phases.
- [ ] Compare a velocity only initialization with an initialization that also
  sets persistent turbulent state before adding such a control.

## TDC

`<prefix>_tdc_face_audit.data` writes active and stored face states, `mlt_vc`,
`gradT`, background `L_conv`, closure luminosities, velocity RHS, and inertia.

Required checks:

- [ ] Check face placement for `mlt_vc`, `Hp`, `gradT`, and `L_conv`.
- [x] Compare the luminosity and velocity rows with `turb:set_TDC_LNA` and the
  nonlinear TDC hydro path.
- [x] Confirm from the hydro source that `v*Uq` and `Eq` are second order on a
  static background and do not add a first-order gas-energy source.
- [x] Keep the nonlinear explicit `mlt_vc_old` operator split out of the static
  LNA eddy-viscosity coefficient.
- [ ] Check turbulent inertia on a static saved model.

The TDC LNA does not differentiate the nonlinear enthalpy flux limiter. The LNA
and nonlinear TDC hydro paths both use the velocity form eddy viscosity
coefficient. They use the same active mixing length. For positive
`harmonic_dissipation_length_beta`, this is the harmonic combination of
`mixing_length_alpha*Hp_hse` and `beta*r`.

GYRE schema 130 does not yet carry the harmonic length. MESA warns and exports
the legacy `mixing_length_alpha` and `Hp_face` values, so a GYRE calculation
from that file uses `Lambda = mixing_length_alpha*Hp_face`.

## MLT And Frozen Flux

Required checks:

- [ ] Run an `mlt_static` saved model and inspect selected modes and work output.
- [ ] Confirm that a `frozen_flux` row structure has no internal TDC `w` rows.
- [ ] Confirm `delta Lconv = 0` from the frozen flux luminosity perturbations.
- [ ] Define an RSP2 frozen flux equation before allowing `frozen_flux` with
  RSP2.

## Boundary Conditions

Required checks:

- [x] Implement a default-off `star_LNA_T_inner` control. The active domain is
  selected by `setup_star_LNA_var_map`, and
  `add_ad_partials_to_matrix` fixes all perturbations below the cut to zero.
- [ ] Check the fixed-core radius and luminosity perturbations selected by
  `star_LNA_T_inner`.
- [ ] Compare envelope kicks for `v_flag` and `u_flag`, including the
  reconstructed face at the core boundary.
- [ ] Check surface momentum and pressure rows for each accepted outer boundary
  condition.
- [ ] Compare `use_RSP_L_eqn_outer_BC` and `RSP2_use_L_eqn_at_surface` with the
  nonlinear branch selection.
- [ ] Check innermost pressure work and face mass for envelope models.
- [ ] Add full star center conditions before applying this operator to full star
  models.

## Solver Cost

The dense path uses DGGEV after algebraic elimination. `star_LNA_num_modes`
limits selected output; it does not reduce the eigensolve dimension. A sparse,
banded, or iterative solver is required before using large zone counts where a
dense all roots solve is too expensive.

The full and reduced scaling passes and each tested full-pencil residual cost
quadratic work in the matrix dimension. The dense elimination and all-root
eigensolve cost cubic work, so the new validation does not change the
asymptotic scaling. It can still add visible wall time through eight reduced
scaling passes and one full residual evaluation per tested candidate. A
convergence stop for equilibration and a row-sparse residual evaluation are
possible local optimizations. The substantial speedup requires a targeted
sparse or iterative solver, or a smaller envelope selected with
`star_LNA_T_inner`.

The first validation target is an envelope or reduced mesh saved model. A
targeted sparse, banded, or iterative solver is required before treating large
MESA/star meshes as a supported workload.

`star_LNA_num_modes` limits selected output and velocity kick indices. It does
not reduce the `DGGEV` matrix dimension.

## Deferred Interfaces

These additions are not required for the first saved model comparisons:

- A standalone executable controlled by an inlist. Such an executable would still need
  to initialize MESA/star, load the model, call star LNA, write output, and
  terminate.
- A selection window based on period or pulsation constant. Raw roots and eigenfunctions
  must first establish that the operator produces the expected acoustic modes.
- Separate decomposed work output and output compatible with RSP. Matched output is
  required before defining another convention.
- An RSP2 initialization control that writes persistent turbulent variables as
  well as velocity. A velocity only comparison must precede this addition.
- Additional source files for solver and output code. The present split is
  sufficient for the saved model audit.

## Rejected Model Options

`check_star_LNA_model` rejects

- rotation;
- RTI;
- mass corrections;
- `other_momentum`, `other_pressure`, and `other_surface_PT` hooks;
- velocity drag;
- `use_compression_outer_BC`;
- the eps_grav energy equation form;
- `constant_L`;
- RSP2 eddy viscosity combined with `u_flag`; and
- nonzero `mstar_dot`.

The LNA supports `u_flag` through the Riemann face reconstruction and cell
momentum equation. It supports the `dPrad/dm` temperature-gradient form with
the hydro opacity floor, face reconstruction, MLT or RSP2 radiative luminosity,
and radiative flux limiter. Perturbed TDC derives `gradT_actual` from this transport
relation before evaluating `Lrad + Lconv - L`. Initial velocity kicks use the
active `u` or `v` hydro variable.
Active RSP2 or MLT turbulent pressure with `u_flag` requires
`star_LNA_perturb_turbulent_pressure = .true.`.

The TDC LNA reports, but does not reject, the nonlinear enthalpy flux limiter,
velocity time centering, and artificial viscosity pressure. The continuous
`u_flag` operator does not apply finite-step velocity time centering.

## Cepheid simple-work correction

The 2026-08-29 `dev_TDC_Cepheid_Pulsation` run selected periods of 22.04 and
9.94 days, while the nonlinear cycle detector measured 11.71 and 12.23 days.
The selected eigenvector residuals were between `7.9d-3` and `4.2d-2`. These
roots are not accurate enough to use for an initial velocity perturbation. A
previous 300-zone run found an 11.826-day fundamental with a residual of
`3.5d-8`.

The work control had been changed from a time-centering-only switch to
`use_P_d_1_div_rho_form_of_work`, but the Star LNA energy row still used the
face-pressure flux for both branches. The nonlinear simple form is

```math
\left(\frac{dW}{dm}\right)_k
= \frac{P_k}{\Delta m_k}
  \left(A_k v_k-A_{k+1}v_{k+1}\right).
```

For a static background, its first-order perturbation is

```math
\delta\left(\frac{dW}{dm}\right)_k
= \frac{P_{k,0}}{\Delta m_k}
  \left[A_{k,0}\,\delta v_k-A_{k+1,0}\,\delta v_{k+1}\right].
```

`dwork_dm_for_star_LNA` now uses this cell-pressure expression when the simple
work control is true. The face-pressure expression remains in the non-simple
branch. This matches the branch used by `hydro_energy:eval_simple_PdV_work` and
the corresponding internal-energy inertia selected by
`mechanical_energy_inertia_for_star_LNA`.

The subsequent equilibrated solve recovered an unstable 11.8309-day raw root,
consistent with the nonlinear period. All otherwise selectable roots still
failed the `1d-3` full-pencil residual limit. The same run exposed a separate
physical mismatch: perturbed TDC used the QHSE spatial gradient despite the
active `dPrad/dm` equation, and its maximum normalized background luminosity
residual was `0.426`. This transport mismatch is corrected in source. The
selected-mode residual check remains the final eigenvector validity test, and
the next run will report the maximum-residual row for the first rejected mode.

Repeating the calculation with `use_dPrad_dm_form_of_T_gradient_eqn = .false.`
left the unstable 11.8313-day fundamental in the raw roots but rejected it from
the selected list. In this branch, the QHSE spatial gradient used interpolated
cell pressure while the TDC radiative coefficient used reconstructed face
pressure. Their relative difference at zone 298 was `2.476d-3`, consistent
with the `2.392d-3` maximum normalized luminosity-row residual. The shared
QHSE helper now uses reconstructed pressure, matching the nonlinear structure
equation and the rest of the TDC face state.

The first post-fix calculation reduced the background luminosity residual to
`9.015d-9`, but the 11.8313-day root had a `1.557d-2` full-pencil residual in
the energy row at `k = 348`. The solver had retained all 350 `w` variables in
the dynamic pencil even though 228 zones impose the algebraic constraint
`delta w = 0`. `partition_star_LNA_indices` now determines dynamic rows from
the presence of a time derivative.

The next calculation confirmed the expected partition: 1172 dynamic
variables, 928 eliminated algebraic variables, and 122 active TDC `w`
equations. It recovered the unstable 11.8313-day fundamental with
`logKE_per_cycle = 8.209d-2`. Its full-pencil residual was `2.796d-2`, so the
partition correction did not remove the remaining numerical defect. A later
model-200 calculation with 118 active TDC zones used 1168 dynamic variables
and reduced the fundamental residual to `4.395d-3` while retaining the same
11.83-day mode. This variation with the background state is consistent with
conditioning or reconstruction error, not incorrect acoustic-mode selection.

Before changing the eigensolver or residual limit, report the reduced-pencil
residual before full eigenvector reconstruction. A small reduced residual with
a larger full residual would localize the loss to the explicit Schur
complement or algebraic reconstruction.

Checks still required:

- [x] Compile MESA and rerun the model-200 Star LNA calculation after the
  `dPrad/dm` and QHSE face-pressure corrections.
- [x] Require selected-mode residuals below `1d-3` before enabling the kick.
- [x] Confirm that the fundamental returns near the measured 12-day period.
- [x] Confirm that `max_abs_luminosity_row_resid` drops from `0.426` to the
  background solve tolerance.
- [x] Repeat with `dPrad/dm` disabled and confirm that the `2.392d-3`
  luminosity-row mismatch is removed.
- [x] Eliminate zero-`w` algebraic constraints from the dynamic pencil.
- [x] Complete a fresh MESA install after the algebraic partition change.
- [x] Confirm 1172 dynamic variables and 928 algebraic variables when 122 TDC
  `w` zones are active.
- [ ] Reduce the full-pencil residual below `1d-3` for the 11.83-day root.
- [ ] Compare reduced-pencil and reconstructed full-pencil residuals for the
  same selected root.
- [ ] Run the full-pencil refinement with a `1d-8` acceptance limit and verify
  the residual, period, growth rate, and kick eigenfunction.
- [ ] Repeat with the LNA delayed until the post-remesh background has relaxed
  for at least one pulsation period.
