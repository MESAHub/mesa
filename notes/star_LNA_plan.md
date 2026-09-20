# star_LNA Implementation Notes

Status: synchronized with `origin/EbF/star_lna` at `622075fbf` on 2026-09-18,
including `origin/main` at `fd396fd73`. The original source audit was performed
on 2026-09-16 at `f6d606939`; the surface-boundary correction is now included.
The existing run record includes a successful 2026-08-29 install after the
TDC transport and algebraic-partition corrections. The model-200 runs recovered
the 11.83-day fundamental but did not meet the full-pencil residual limit.
Validation of Newton refinement and the later August 31/September 1 changes
remains open in `star_LNA_readiness_checklist.md`. No MESA build or model run
was performed for this documentation update.

The manuscript is `star_LNA_manuscript.tex`; its PDF is the continuous
derivation of the current local prescription. Dated audit and development
entries below preserve the earlier implementation history. The 2026-09-16
entry records the source changes and the surface-boundary fix that was then
only on the remote branch. The 2026-09-18 entries record the equation layout,
the restored energy-form and velocity-grid distinctions, and upstream
synchronization with local work preserved.

These notes replace the repository root `star_LNA_plan.md`. They document the
implementation that is in the tree now, not the original proposal. The separate
rewrite plan lives in `notes/star_LNA_compression_plan.md`.

## Scope

`star_LNA` performs radial linear nonadiabatic analysis of a static MESA/star
model. It uses MESA/star variables and discretizations. RSP LINA defines the
output layout, growth rate convention, velocity kicks, and work terms used for
comparison.

Implemented paths:

- A face velocity LNA operator. Hydro `v_flag` may be inactive because the LNA
  introduces its own linear velocity perturbation.
- A cell velocity operator for `u_flag`, using the MESA Riemann face
  reconstruction and cell momentum equation without finite-step time
  centering.
- Dense generalized eigenproblem solves with algebraic elimination and local
  full-pencil Newton refinement of candidate eigenpairs.
- A temperature-selected envelope domain through `star_LNA_T_inner`.
- An explicit problem object containing the variable map, equation registry,
  and dense matrices.
- Base variables `lnd`, `lnR`, active radial velocity (`v` or `u`), `lnT`,
  and `L`.
- RSP2 variables `w` and `Hp` when `RSP2_flag` is active.
- An internal TDC `w` variable when MLT/TDC is active and convection is
  perturbed.
- Perturbed radiative luminosity and convective luminosity.
- Turbulent pressure and turbulent energy terms for RSP2.
- Eddy viscosity acceleration for RSP2 and TDC.
- RSP period/growth sorting, mode files, work files, and optional velocity
  kicks from as many as three selected modes.
- Row structure with named equations, RSP2 term audit, and TDC face state audit output.
- A startup or model number analysis path that can stop before evolution after a
  successful solve.

Rejected model options:

- Rotation and corrections for rotation.
- RTI source, diffusion, and energy terms.
- Gravitational/baryonic mass corrections.
- User `other_*` hooks, including other momentum, pressure, and surface PT.
- Velocity drag.
- `use_compression_outer_BC`.
- The `eps_grav` energy equation form.
- `constant_L`.
- Nonzero `mstar_dot`.
- RSP2 eddy viscosity with `u_flag`. The nonlinear Riemann path still places
  this term in the face reconstruction rather than the cell momentum source.

`check_star_LNA_model` rejects these options before matrix assembly.
Both the logarithmic-pressure and `dPrad/dm` temperature-gradient forms are
implemented. Static MLT and RSP2 use the selected temperature-gradient row
directly. Perturbed TDC obtains its spatial gradient from the selected equation
and retains the coupled `Lrad + Lconv - L` closure. Saved model validation of
the static velocity MLT, TDC, RSP2, and `u_flag` operators remains.
The required runs are listed in `star_LNA_readiness_checklist.md`.

For `u_flag`, active RSP2 or MLT turbulent pressure requires
`star_LNA_perturb_turbulent_pressure = .true.` because the Riemann face state
includes those pressure terms.

The setup report identifies accepted controls that the LNA does not
differentiate. These include the TDC enthalpy flux limiter, velocity time
centering, and artificial viscosity pressure. TDC eddy viscosity uses cell stresses for face
velocity and face stresses for cell velocity, each formed from the current
static convection state.

## Source Map

Primary implementation:

- `star/private/star_LNA.f90`, public entry point, top level flow, and
  main equation assembly order.
- `star/private/star_LNA_support.f90`, variable map, row implementations,
  solver, output, work/kick, and audit helpers.
- `star/private/star_LNA_turbulence_closures.f90`, pressure, convection,
  turbulent energy, and eddy-viscosity closures.
- `turb/public/turb.f90`, public wrapper `set_TDC_LNA`

Integration and controls:

- `star/Makefile`
- `star/job/run_star_support.f90`, `maybe_do_star_LNA`
- `star/public/star_lib.f90`
- `star/defaults/controls_dev.defaults`
- `star_data/private/star_controls_dev.inc`
- `star/private/ctrls_io.f90`

RSP reference implementation:

- `star/private/rsp_lina.f90`
- `star/private/rsp.f90`
- `star/private/rsp_build.f90`

Local papers and extracted notes are kept outside the branch under
`Classical_Pulsations/star_lna/notes/references/`:

- `smolec_phd_thesis.pdf` and its text extraction.
- `smolec_moskalik_2008_convective_hydrocodes.pdf` and its text extraction.
- `farag_etal_2026_arxiv_2603.15766_source/`.

The Smolec text files are Poppler `pdftotext` extractions. Use them to locate
material. Verify equations against the PDFs and MESA source.

## Reference: What RSP LINA Does

`rsp_lina.f90` is self contained and linearized by hand. Its state vector is

```text
X = {dR_1, dU_1, dT_1, dw_1, ..., dR_N, dU_N, dT_N, dw_N}
```

and it builds an ordinary eigenproblem:

```text
LLL * X = SIGMA * X
```

with four rows per zone:

```text
dR/dt = U
dU/dt = -4*pi*R^2*dP/dm - G*M/R^2
cv*dT/dt = -(P + (de/dV)_T)*dV/dt - dL/dm
de_t/dt = source - dissipation - dLt/dm - P_t*dV/dt
```

Important RSP conventions to preserve in star LNA:

- Normal modes use `exp(sigma*t)`.
- `sigma = sigma_re + i*sigma_im`.
- `period = 2*pi/sigma_im`.
- The printed RSP LINA growth is the log kinetic energy growth per period.
  Star LNA labels this as `logKE/cyc` in terminal output and
  `logKE_per_cycle` in data files to avoid confusion with GYRE's different
  `eta` convention:

```text
logKE_per_cycle = 4*pi*sigma_re/sigma_im
KE_fractional_growth = exp(logKE_per_cycle) - 1
GREKM = 2*tanh(logKE_per_cycle/2)
amplitude_fractional_growth = exp(logKE_per_cycle/2) - 1
```

- RSP sorts by increasing positive `sigma_im`.
- RSP skips very strongly damped low frequency entries using practical growth
  filters:

```text
first accepted mode: growth > -3
later accepted modes: growth > -5
```

- Star LNA also applies `star_LNA_min_mode_frequency_uHz` before the growth
  filters when selecting the summary modes. The raw finite positive frequency
  list is still written without this filter. The default `1d-3` microHz cutoff
  removes very slow thermal/quasistatic roots from the selected output while
  preserving nearly neutral acoustic modes with `eta` as small as `1d-8`.

- RSP writes period/growth, eigenfunction, and work integral files.
- RSP kicks nonlinear runs with a real velocity profile built from the radial
  displacement eigenfunction and normalized at the surface.

Smolec and Moskalik 2008 describe the same convention, using `eta` for
`4*pi*Re(sigma)/omega`, where `omega = Im(sigma)`. They also note that the
added turbulent energy equation creates a separate branch of strongly damped
turbulent modes. Account for this branch when classifying MESA/star roots. RSP
defines comparison behavior, not the MESA/star row discretization.

## star_LNA Matrix Convention

The star LNA code builds a generalized eigenproblem:

```text
A * dx = sigma * B * dx
```

`sigma` has units `1/s`. A continuous equation

```text
dq/dt = F(q)
```

is assembled as

```text
dF = sigma * dq
```

so the partials of `F` go in `A` and the partials of the inertial quantity `q`
go in `B`.

An algebraic equation

```text
G(q) = 0
```

is assembled only into `A`.

The dense solver path scales the full rows and columns, partitions and
eliminates the algebraic variables, and equilibrates the reduced generalized
pencil. It divides the reduced `B` matrix by the model dynamical timescale so
`DGGEV` solves for `lambda = sigma*t_dyn`, then reconstructs and unscales the
full eigenvectors. Candidates above `star_LNA_max_eigenvector_residual` enter
full-pencil Newton refinement if their componentwise residual is at most
`1d-1` or their balanced normwise residual is at most `1d-6`. The refined pair is
retained only if it reduces the residual against the original unscaled
equations. The acceptance limit and mode filters are then applied. The
equations and banded solves are given in the full-pencil refinement section.

There is no fixed star_LNA matrix size cutoff. The dense path is limited by
available memory and by the cost of the dense LAPACK solve.

Full-pencil scaling and eight reduced-pencil equilibration passes add
quadratic work. Candidate residual checks and Newton refinement use the full
matrix bandwidth. Dense elimination and `DGGEV` remain cubic and dominate as
the active zone count grows. A targeted sparse or iterative solver is needed
to avoid computing the full reduced spectrum.

`star_LNA_num_modes` does not make the dense solve cheaper. `DGGEV` computes the
full reduced eigensystem first, and `star_LNA_num_modes` only limits how many
selected modes are written or used for the kick. This is the same scaling issue
RSP LINA would have at the same zone count; RSP is usually cheaper because its
models are typically much smaller.

Algebraic elimination is:

```text
Aaa*x_alg + Aad*x_dyn = 0
x_alg = -inv(Aaa)*Aad*x_dyn

(Add - Ada*inv(Aaa)*Aad)*x_dyn =
   sigma*(Bdd - Bda*inv(Aaa)*Aad)*x_dyn
```

In code:

- `build_reduced_star_LNA_problem` partitions the variables.
- `DGESVX` equilibrates and refines the solve for `inv(Aaa)*Aad`.
- `DGGEV` solves the reduced problem.
- `reconstruct_star_LNA_eigenvectors` applies the minus sign when restoring
  algebraic components.

## Current Top Level Flow

`do_star_LNA` currently does:

```fortran
call set_vars_if_needed(s, 0d0, 'star_LNA', ierr)
call check_star_LNA_model(s, ierr)
call setup_star_LNA_problem(s, problem, ierr)
call report_star_LNA_setup(s, problem%map)
call assemble_star_LNA_equations(s, problem, ierr)
call write_star_LNA_matrix_summary(s, problem, ierr)
call solve_dense_star_LNA(s, problem%map, problem%mtx, ierr)
call free_star_LNA_problem(problem)
```

`do_star_LNA` and `assemble_star_LNA_equations` live in
`star/private/star_LNA.f90`. The main file therefore shows both the public solve
flow and the main equation ordering. The detailed row stencils, closures,
solver/output logic, and diagnostics live in `star/private/star_LNA_support.f90`.

`assemble_star_LNA_equations` calls the row assemblers in this order:

```fortran
call assemble_density_rows
call assemble_radius_rows
call assemble_momentum_rows
call assemble_energy_rows
call assemble_luminosity_rows
call assemble_rsp2_turbulent_rows
call assemble_tdc_turbulent_rows
```

This flow is correct enough to keep, but the implementation below it is too
large and too spread out. The compression plan proposes keeping this main
shape while replacing the internal organization.

## Variable Map

Base variables per zone:

```text
lnd  = delta ln rho
lnR  = delta ln R
v    = delta velocity
lnT  = delta ln T
L    = delta luminosity
```

RSP2 adds:

```text
w    = delta turbulent velocity / turbulent variable
Hp   = delta pressure scale height
```

TDC, when active, adds:

```text
w    = internal LNA TDC turbulent velocity variable
```

The TDC `w` is not a new hydro solver state. It exists only inside the LNA
eigenproblem so the TDC velocity relation can contribute a B matrix inertia
instead of being forced into a purely algebraic luminosity perturbation.

Dynamic variables:

```text
lnR, v, u, lnT, w
```

Algebraic variables:

```text
lnd, L, Hp
```

The active radial velocity variable is `v` for the face momentum equation and
`u` for `u_flag`. Automatic differentiation velocity components map to the
active variable. The radius and pressure-work rows use reconstructed face
velocity in the `u_flag` path. The RSP-compatible velocity kick remains limited
to `v_flag` because it writes `v`, `xh(i_v,:)`, and `v_start`.

## Controls

Current dev controls:

```fortran
star_LNA_flag
star_LNA_model_number
star_LNA_stop_after_run
star_LNA_num_modes
star_LNA_min_mode_frequency_uHz
star_LNA_min_first_mode_eta
star_LNA_min_mode_eta
star_LNA_max_abs_mode_eta
star_LNA_max_eigenvector_residual
star_LNA_output_directory
star_LNA_output_file_prefix
star_LNA_write_matrix_summary
star_LNA_write_period_growth
star_LNA_write_eigenfunctions
star_LNA_write_work_integrals
star_LNA_set_initial_velocity
star_LNA_kick_vsurf_km_per_sec
star_LNA_mode_for_period
star_LNA_kick_mode_1
star_LNA_kick_mode_2
star_LNA_kick_mode_3
star_LNA_kick_fraction_1
star_LNA_kick_fraction_2
star_LNA_kick_fraction_3
star_LNA_solver
star_LNA_include_tdc
star_LNA_include_rsp2
star_LNA_perturb_convective_flux
star_LNA_perturb_turbulent_pressure
star_LNA_perturb_eddy_viscosity
star_LNA_perturb_turbulent_energy
star_LNA_convection_treatment
```

Defaults are in `star/defaults/controls_dev.defaults`:

```fortran
star_LNA_flag = .false.
star_LNA_model_number = -1
star_LNA_stop_after_run = .false.
star_LNA_num_modes = 3
star_LNA_min_mode_frequency_uHz = 1d-3
star_LNA_min_first_mode_eta = -3d0
star_LNA_min_mode_eta = -5d0
star_LNA_max_abs_mode_eta = 10d0
star_LNA_max_eigenvector_residual = 1d-3
star_LNA_output_directory = ''
star_LNA_output_file_prefix = ''
star_LNA_write_matrix_summary = .true.
star_LNA_write_period_growth = .true.
star_LNA_write_eigenfunctions = .true.
star_LNA_write_work_integrals = .true.
star_LNA_set_initial_velocity = .false.
star_LNA_kick_vsurf_km_per_sec = 1d0
star_LNA_mode_for_period = 0
star_LNA_kick_mode_1 = 1
star_LNA_kick_mode_2 = 0
star_LNA_kick_mode_3 = 0
star_LNA_kick_fraction_1 = 1d0
star_LNA_kick_fraction_2 = 0d0
star_LNA_kick_fraction_3 = 0d0
star_LNA_solver = 'dense'
star_LNA_include_tdc = .true.
star_LNA_include_rsp2 = .true.
star_LNA_perturb_convective_flux = .true.
star_LNA_perturb_turbulent_pressure = .true.
star_LNA_perturb_eddy_viscosity = .true.
star_LNA_perturb_turbulent_energy = .true.
star_LNA_convection_treatment = 'perturbed'
```

`star_LNA_model_number < 0` runs at startup. A nonnegative value runs when
`s%model_number` matches in `maybe_do_star_LNA`.

`star_LNA_stop_after_run = .true.` stops the normal star driver before
evolution after a successful analysis. It is intended for loaded model
diagnostics with `star_LNA_set_initial_velocity = .false.`. Leave it false for
velocity kick runs that continue into nonlinear evolution. The binary driver
does not implement this stop workflow and rejects a startup request rather than
silently continuing binary evolution.

`star_LNA_solver = 'dense'` is the only implemented solver. It computes the full
reduced eigensystem before applying `star_LNA_num_modes`, so lowering
`star_LNA_num_modes` reduces output and kick work but not the LAPACK solve cost.

`star_LNA_min_mode_frequency_uHz` is a selection filter only. It does not alter
the matrix or dense solve. Set it to `0d0` to disable the frequency floor, or
raise it for Solar pressure mode output.

`star_LNA_min_first_mode_eta` and `star_LNA_min_mode_eta` expose the old
fixed in source RSP eta filters. `star_LNA_max_abs_mode_eta` is an optional
symmetric guard against pathological roots with huge growth or decay per
period; set it to `0d0` to disable.

`star_LNA_max_eigenvector_residual` is the largest componentwise backward
error accepted against the original full pencil. It is an upper limit, not a
floor. It must be positive, so the eigenvector validity check cannot be
disabled. The default is `1d-3`; `1d-1` is useful only for inspecting rejected
roots because it permits a ten percent row defect.

### Convection Treatment Semantics

The LNA controls intentionally separate a local "fixed convective luminosity"
switch from the global row choice:

- `star_LNA_perturb_convective_flux = .false.` means the active branch should
  hold its convective luminosity expression fixed where that branch supplies one.
  It is a derivative switch for `Lconv`, not a request to select another
  luminosity equation.
- `star_LNA_convection_treatment = 'perturbed'` uses the active dynamic
  convection path. For TDC this adds the internal LNA `w` variable and uses
  the local TDC `Lrad+Lconv-L` closure plus the TDC velocity row.
- `star_LNA_convection_treatment = 'frozen'` is the older broad diagnostic mode.
  For TDC models without RSP2, it disables `tdc_lna_active`, removes the internal TDC
  `w` variable, and lets the luminosity slot fall back to the surface temperature
  row at `k = 1` and the static temperature gradient row for `k > 1`. This is
  not the same as a GYRE frozen convective flux calculation.
- `star_LNA_convection_treatment = 'frozen_flux'` selects the frozen convective
  flux row. It keeps the background convective luminosity in
  the equilibrium relation but sets its perturbation to zero:

```text
L_resid = Lrad_ad + Lconv0 - L_ad
delta Lconv = 0
```

`frozen_flux` does not add the internal TDC `w` variable. RSP2 rejects this
value because its frozen flux equation is not defined.

- `star_LNA_convection_treatment = 'mlt_static'` explicitly selects the
  static MLT temperature gradient row for models without TDC or RSP2. It adds
  no convective velocity variable and is rejected for TDC and RSP2 backgrounds.

`star_LNA_set_initial_velocity` applies a kick from the selected complex
velocity eigenfunction(s). Each mode is phase-aligned and normalized by its
surface velocity before the requested fractions are combined. The result is
scaled to the requested surface velocity. The LNA solve can run without a
hydro velocity variable. The kick requires `v_flag` or `u_flag` and writes the
active hydro variable and its start state. The `u_flag` path also rebuilds the
Riemann face velocity and initializes `u_face_start`. The kick is rejected when
`use_fixed_vsurf_outer_BC` is active, because that surface BC constrains the
LNA surface velocity and makes `star_LNA_kick_vsurf_km_per_sec` undefined.
After a kick, star_LNA reports the maximum absolute velocity and maximum
`|velocity|/cs`; a large interior/surface amplitude ratio points to a bad kick
eigenfunction or mode choice rather than a surface BC conversion issue.

## Base Equation Rows

The rows below use the current code's variable choices and sign convention:

```text
A*dx = sigma*B*dx
```

### Density / Volume Closure

Code:

- `assemble_density_rows`
- `cell_volume_for_star_LNA`

The row enforces the Lagrangian cell mass/volume closure:

```text
rho_k * dV_k = dm_k
```

The logarithmic residual is:

```text
G_rho,k = ln(rho_k) + ln(Delta V_k) - const
```

The constant is irrelevant for linearization, so the row is:

```text
d ln rho_k + d ln Delta V_k = 0
```

with

```text
Delta V_k = (4*pi/3)*(r_k^3 - r_{k+1}^3)
```

For the innermost envelope cell, `r_{k+1}` is the fixed center radius through
`wrap_r_p1`.

### Radius / Kinematics

Code:

- `assemble_radius_rows`
- `force_zero_velocity_for_star_LNA`

For normal face velocity LNA cells:

```text
d ln R_k/dt = v_k/r_k
```

Linearized:

```text
d(v_k/r_k) = sigma*d ln R_k
```

So `d(v/r)` is inserted into `A`, and `d lnR_k` is inserted into `B`.

If the hydro setup force zeros the velocity, the LNA row is:

```text
d v_k = 0
```

### Momentum

Code:

- `assemble_momentum_rows`
- `momentum_rhs_for_star_LNA`
- `surface_velocity_rhs_for_star_LNA`
- `dPtot_face_for_star_LNA`
- `d_mlt_Pturb_face_for_star_LNA`
- `Uq_face_for_star_LNA`
- `star_LNA_dm_face`

For interior faces, the current row is:

```text
dv_k/dt = grav_k + Uq_k
          - (Delta Ptot_k + Delta Pmlt_turb_k)/(dm_face_k/A_k)
```

with

```text
grav_k = -G*m_k/r_k^2
A_k    = 4*pi*r_k^2
```

Linearized:

```text
d[grav + Uq - A_k*(Delta Ptot + Delta Pmlt_turb)/dm_face]
   = sigma*d v_k
```

`Ptot` currently means EOS pressure plus RSP2 turbulent pressure when RSP2
turbulent pressure perturbations are active. MLT/TDC turbulent pressure is added
through `d_mlt_Pturb_face_for_star_LNA`, because that pressure is face centered
and uses MLT/TDC convective velocity rather than RSP2 cell `w`.

TDC MLT turbulent pressure uses the internal LNA `w` when TDC LNA is active:

```text
P_mlt_turb = mlt_Pturb_factor * rho_face * conv_vel^2/3
```

The eddy viscosity acceleration `Uq` now uses a local static `Chi` helper for
RSP2, TDC, and the fallback path. It does not reuse cached hydro `Chi_ad` or
normal TDC `compute_Chi_cell`, because those helpers can carry nonlinear
background velocity gradient terms. In the static LNA, only the coefficient
times `delta(v/r)` is first order.

For the TDC closure, write the cell stress as

```text
Chi = C*w*Delta(v/r).
```

Since `Delta(v/r) = 0` in the static background,

```text
delta Chi = C0*w0*delta[Delta(v/r)].
```

The `TDC_alpha_M_use_explicit_mlt_vc_in_momentum_equation` control replaces
the current `mlt_vc` with the start-of-step `mlt_vc_old` only as an operator
split in the nonlinear momentum solve. It does not define a different
continuous linear operator. The TDC eddy-viscosity coefficient in `star_LNA`
therefore uses the current equilibrium `s%mlt_vc` whether or not the nonlinear
explicit control is enabled.

The face mass is selected by `star_LNA_dm_face`, which was changed to follow the
MESA hydro momentum face convention rather than assuming `s%dm_bar(k)` is always
the correct denominator.

The outer row can be a fixed surface velocity row, a pressure boundary row, or a
surface momentum row depending on active controls. `use_RSP_L_eqn_outer_BC` is
handled in the luminosity row, not the momentum row.

### Energy

Code:

- `assemble_energy_rows`
- `energy_rhs_for_star_LNA`
- `dL_dm_for_star_LNA`
- `energy_sources_for_star_LNA`
- `dwork_dm_for_star_LNA`
- `energy_inertia_for_star_LNA`
- `mechanical_energy_inertia_for_star_LNA`
- `turbulent_energy_inertia_for_star_LNA`

The current gas energy row is:

```text
d e_eff,k/dt = -dL_k/dm + sources_k - dwork_k/dm
```

Linearized:

```text
d[-dL/dm + sources - dwork/dm] = sigma*d e_eff
```

When `use_P_d_1_div_rho_form_of_work` is true, the static work perturbation is

```text
delta(dwork_k/dm) = P_k*[delta(A*v)_k - delta(A*v)_{k+1}]/dm_k.
```

The non-simple total-energy branch retains the face-pressure flux difference.

The luminosity divergence is:

```text
dL/dm = (L_k - L_{k+1})/dm_k
```

The non-simple pressure work stencil uses static face pressure multiplying the
linearized velocity divergence:

```text
dwork/dm =
  4*pi*(P_out*r_k^2*d v_k - P_in*r_{k+1}^2*d v_{k+1})/dm_k
```

For the innermost envelope cell:

```text
dwork/dm = 4*pi*P_out*r_k^2*d v_k/dm_k
```

This is RSP in the sense that pressure perturbations do not multiply a
background compression rate in the static limit. It is also compatible with
MESA's continuous first order linearization when the background velocity is zero.
If we ever choose to linearize a nonstatic background, this row must be revisited
because `dP*v0` and `P*dA*v0` terms would then be first order.

The inertial energy is:

```text
d e_eff =
  (dE/dRho)*rho*d ln rho + Cv*T*d lnT
  + optional d e_turb
  + optional d e_kinetic + optional d e_potential
```

The mechanical energy pieces are included only when the active MESA form is the
`dedt` total-energy branch and `use_P_d_1_div_rho_form_of_work` is false. In the
static limit `d e_kinetic` vanishes, but `d e_potential` is first order in the
radius perturbation.

RSP2 turbulent inertia:

```text
e_turb = etrb = w^2
```

TDC turbulent inertia, when
`TDC_include_eturb_in_energy_equation` and `star_LNA_perturb_turbulent_energy`
are active:

```text
e_turb = 0.75*(mlt_vc_k^2 + mlt_vc_{k+1}^2)
```

using the internal TDC `w` path when TDC LNA is active.

Current limitation:

```text
turbulent_viscous_heating_for_star_LNA = 0
```

So `Eq` heating is not inserted as a first order source in the gas energy row.
This matches a static background expectation for quadratic velocity gradient
heating, but it is still an audit item because MESA's TDC and RSP2 bookkeeping
also use `Eq` in related turbulent equations and diagnostics.

### Luminosity / Temperature Row

Code:

- `assemble_luminosity_rows`
- `rsp_lsurf_resid_for_star_LNA`
- `tdc_luminosity_resid_for_star_LNA`
- `surface_temperature_resid_for_star_LNA`
- `temperature_gradient_resid_for_star_LNA`

The row selection is:

```fortran
if use_rsp_lsurf_row_for_star_LNA(s, k):
   L - L_surf = 0
else if tdc_lna_active(s) and k > 1:
   Lrad + Lconv - L = 0
else if frozen_flux_lna_active(s) and k > 1:
   Lrad + Lconv0 - L = 0
else if k == 1:
   surface T boundary = 0
else:
   temperature gradient residual = 0
```

RSP2 uses the selected temperature-gradient residual. Its local
`hydro_rsp2:compute_RSP2_gradT` closure enforces `Lr + Lc + Lt = L` before
that row is assembled. There is no separate RSP2 luminosity-sum row.

For TDC:

```text
d(Lrad_k + Lconv_k - L_k) = 0
```

For frozen convective flux:

```text
d(Lrad_k - L_k) = 0, with Lconv0 held fixed
```

For static MLT or RSP2:

```text
d[gradT_from_temperature_difference - gradT_ad] = 0
```

The current TDC luminosity row applies only for `k > 1`. Surface behavior still
uses the surface T boundary unless the RSP surface luminosity row is selected.

## RSP2 Rows

Code:

- `assemble_rsp2_turbulent_rows`
- `assemble_rsp2_w_row`
- `rsp2_turbulent_energy_rhs_for_star_LNA`
- `rsp2_turbulent_pressure_work_dm_for_star_LNA`
- `rsp2_dLt_dm_for_star_LNA`
- `rsp2_turbulent_energy_inertia_for_star_LNA`
- `assemble_rsp2_Hp_row`

RSP2 adds a cell turbulent-energy row and an independent face flux row.
Its luminosity and source terms use the same signed `s%PII_ad` relation:
`PII = x_ALFAS*(Lambda_face/Hp_face)*Cp_face*Y_face`.
`compute_RSP2_gradT` sets `gradT = gradL + Y_face` directly.
The local damping in `rsp2_damping_for_star_LNA` is proportional to `w^3`;
it contains no `w_min^3` subtraction.

The turbulent energy row is:

```text
d etrb_k/dt = COUPL_k - dLt_k/dm - Ptrb_k*dVdt_k/dm
```

Linearized:

```text
d[COUPL - dLt/dm - Ptrb0*dVdt/dm] = sigma*d etrb
```

The star LNA matrix variable is `w` because that is the native MESA RSP2
variable, but the B side of this row uses `wrap_etrb_00 = w^2`. Therefore the
B matrix derivative is `delta etrb = 2*w0*delta w`, which is equivalent to
RSP LINA's use of `Et` as the turbulent energy variable.

`Ptrb0` is held as the static pressure coefficient in the compression work term
because `dPtrb*dVdt0` vanishes for a static background.
The nonlinear MESA turbulent energy residual has an `Eq` source term, but it is
quadratic in the velocity gradient perturbation and is omitted from the static
first order LNA, matching the RSP LINA treatment.

The turbulent flux divergence is:

```text
dLt/dm = (Lt_k - Lt_{k+1})/dm_k
```

The scale height row is algebraic:

```text
d(Hp_expected - Hp) = 0
```

Cells forced nonturbulent by RSP2 controls use:

```text
d w = 0
```

## TDC Perturbed Convection

Code:

- `turb/public/turb.f90:set_TDC_LNA`
- `tdc_lna_active`
- `tdc_luminosity_resid_for_star_LNA`
- `tdc_luminosity_terms_for_star_LNA`
- `tdc_relation_for_star_LNA`
- `assemble_tdc_turbulent_rows`
- `assemble_tdc_w_row`
- `tdc_A_for_star_LNA`
- `tdc_conv_vel_for_star_LNA`
- `tdc_Eq_div_w_for_star_LNA`

TDC does not normally carry an independent star solver variable for turbulent
velocity in the same way RSP2 does. The current LNA adds an internal `w` variable
so the TDC turbulent velocity relation can be represented as a dynamic row.

The internal variable is:

```text
A = conv_vel/sqrt(2/3)
conv_vel = sqrt(2/3)*A
```

The TDC wrapper returns:

```text
luminosity_resid = L_rad + L_conv - L_total
velocity_rhs
velocity_inertia
L_rad
L_conv
```

The wrapper computes:

```text
L0 = (16*pi*a*c/3)*G*m*T^4/(P*kappa)
gradr = L_total/L0
Lambda0 = alpha_MLT*Hp_for_mlt
Lambda = Lambda0                                      if beta_harm <= 0
Lambda = 1/(1/Lambda0 + 1/(beta_harm*r))             if beta_harm > 0
alpha_eff = alpha_MLT                                 if beta_harm <= 0
alpha_eff = Lambda/Hp_hse                             if beta_harm > 0
```

For harmonic dissipation, `Hp_hse = P/(rho*g)`. Otherwise `Hp_for_mlt` is the
active reconstructed or stored face scale height. The implementation calls
`star_utils:get_mlt_mixing_length` for `Lambda`.

It calls `set_MLT('Cox', ...)` for the local MLT quantities, then uses:

```text
Y = gradT_actual - gradL

if Y > 0 and include_mlt_corr_to_TDC:
   Y_env = Y*Gamma/(1 + Gamma)
else:
   Y_env = Y
```

The optional enthalpy flux limiter is ignored by star LNA.  The LNA target is
the unsaturated local TDC relation, so the source is always:

```text
S0 = TDC_alpha_S*x_ALFAS*alpha_eff*Cp*T/Hp_for_mlt*grada*Y_env
     + Eq_div_w
```

The damping terms are:

```text
D0  = TDC_alpha_D*x_CEDE/Lambda
DR0 = 4*sigma_SB*(TDC_alpha_R*x_GAMMAR/Lambda)^2
      *T^3/(rho^2*Cp*kappa)
```

The TDC velocity relation is assembled as:

```text
velocity_rhs = S0 - A*DR0 - A^2*D0
```

and the B matrix inertia is:

```text
velocity_inertia = 2*A
```

If `TDC_alpha_Pt` is nonzero, the wrapper adds the turbulent pressure
compression contribution to the inertia:

```text
velocity_inertia =
   2*A + A0*TDC_alpha_Pt*(2/3)*rho0/rho
```

The implementation intentionally uses `A%val` and `rho%val` in the extra
coefficient. This places the static coefficient multiplying
`sigma*d(1/rho_face)` on the generalized eigenproblem B side.

The luminosities are:

```text
L_rad  = L0*gradT_actual
L_conv = TDC_alpha_C*alpha*alpha_c*rho*T*Cp*4*pi*r^2*A*Y_env
```

The TDC `w` row is:

```text
d velocity_rhs = sigma*d velocity_inertia
```

Cells without convection and cells forced to zero use:

```text
d w = 0
```

The TDC helper is public from `turb` because it reuses the same local TDC/MLT
constants and branch logic. It does not change normal TDC operation and does not
use a finite timestep.

After the rebase onto face reconstruction, `tdc_relation_for_star_LNA` obtains
`T`, `rho`, `P`, `e`, `Cp`, `chiRho`, `chiT`, `grada`, opacity, `Hp`, and
`gradr` from `get_reconstructed_face_state_ad`. The helper returns reconstructed
AD face quantities when `use_face_reconstruction = .true.` and the stored face
quantities otherwise. The TDC radiative luminosity, turbulent pressure
background density, and TDC eddy viscosity `Hp` stencil follow the same choice.

## Pressure, Work, and Luminosity Diagnostics

Eigenfunction and work output is written in `write_star_LNA_eigenfunctions`.
All star LNA diagnostic files are opened under:

```text
<star_LNA_output_directory>/<star_LNA_output_file_prefix>_*.data
```

If `star_LNA_output_directory = ''`, star LNA uses `log_directory`; if both are
empty, it falls back to `LOGS`.  The chosen directory is created if needed with
the same `folder_exists`/`mkdir` pattern used by normal MESA log output.

The path is relative to the run directory unless an absolute directory is given.
For example, the `dev_TDC_Cepheid_6M` pulse inlist sets
`log_directory = 'LOGS_pulsation'`, leaves `star_LNA_output_directory = ''`, and
sets `star_LNA_output_file_prefix = 'star_LNA'`, so the files are written as
`LOGS_pulsation/star_LNA_*.data`, not directly in the model directory.  Set
`star_LNA_output_directory = '.'` to write directly in the run directory.

If `star_LNA_output_file_prefix = ''`, star LNA uses `star_LNA`.  The prefix is
a file stem, not a path.  Leading/trailing slashes in the prefix are ignored,
and internal slashes are converted to underscores so `star_LNA_output_directory`
is the only control that determines directories.

Files are opened with replace semantics.  Multiple LNA calls should use distinct
file prefixes when their output needs to be kept side by side.

Files:

```text
<prefix>_raw_modes.data
<prefix>_period_growth.data
<prefix>_eigenfunction_<mode>.data
<prefix>_work_<mode>.data
<prefix>_matrix_summary.data
<prefix>_row_structure.data
<prefix>_rsp2_term_audit.data
<prefix>_tdc_face_audit.data
```

`<prefix>_matrix_summary.data` and `<prefix>_row_structure.data` are written
after matrix assembly. The matrix summary includes equation counts, and row
structure names the equation in every row. RSP2 and active TDC cases also write
their applicable term or face audit with matrix summary output. Raw modes,
period/growth, eigenfunction, and work files
are written only after the dense eigensolve succeeds.  Eigenfunction and work
files are written once per selected mode.

Eigenfunction files include:

```text
d lnd
d lnR
d v
d lnT
d L
d w      if present
d Hp     if present
d Lr
d Lc
d Lt
abs/phase dlnR
abs/phase dlnT
abs/phase dL/L0
abs/phase w
```

Work files contain these initial components:

```text
pressure_work
turb_pressure_work
eddy_visc_work
rad_lum_work
conv_lum_work
turb_lum_work
total_work
```

The ordinary and turbulent pressure pieces are intended to follow the Smolec/RSP
form:

```text
W_P = -pi * Im{conjg(dP) * dV}
```

The eddy viscosity diagnostic follows the same static background coefficient used
by the momentum row. The luminosity work terms are currently phase diagnostics
from luminosity divergence perturbations; their sign and normalization still need
to be validated against RSP `LINA_work*.data`.

`star_LNA_pressure_perturbations` now uses the same EOS and RSP2 turbulent
pressure component closure as the operator. TDC MLT turbulent pressure remains
a separate face centered momentum and static work term; a dedicated TDC
turbulent pressure work component remains open.

## Mode Output and Kicks

The finite generalized eigenvalue is:

```text
sigma = (alphar + i*alphai)/beta
```

The printed period and RSP growth are:

```text
period_days = 2*pi/abs(sigma_imag)/86400
pulsation_constant_Q_days = period_days*sqrt((M/Msun)*(Rsun/R)^3)
W_rad_per_sec = abs(sigma_imag)
logKE_per_cycle = 4*pi*sigma_real/sigma_imag
```

The code writes all finite positive frequency raw roots sorted by increasing
`sigma_imag` to `<prefix>_raw_modes.data`.

The selected list printed to the terminal and written to
`<prefix>_period_growth.data` uses the RSP positive frequency and
logarithmic kinetic energy growth filters. It writes `period_days`,
`pulsation_constant_Q_days`, `W_rad_per_sec`, and exact derived comparison
columns:
`ke_fractional_growth_per_period = exp(logKE_per_cycle)-1`,
`grekm_growth_per_period = 2*tanh(logKE_per_cycle/2)`, and
`amplitude_fractional_growth_per_period = exp(logKE_per_cycle/2)-1`. There is
no selector based on Q and no documented period window in the current
implementation.

When matrix summary output is enabled, star LNA also writes
`<prefix>_row_structure.data`. This is a static audit aid: for every full matrix
row it records the row variable, A/B nonzero counts, largest absolute A/B entry,
and the dominant A/B variable/zone. It does not affect the solve.

The kick controls start at one in the selected mode list:

```text
star_LNA_kick_mode_1 = 1
star_LNA_kick_mode_2 = 0
star_LNA_kick_mode_3 = 0
```

`star_LNA_kick_mode_1` must be positive. A zero value disables the second and
third kick components. Active kick modes are limited by `star_LNA_num_modes`,
not by an RSP fixed 15-mode array.

`star_LNA_mode_for_period` starts at zero, matching the existing RSP control.

The kick uses the computed complex velocity eigenfunction:

- loads the complex eigenvector,
- phase-aligns each requested mode to its complex surface velocity,
- normalizes each velocity eigenfunction to unit surface velocity,
- combines up to three modes using normalized fractions,
- scales the final velocity so `v(1)` equals
  `1d5*star_LNA_kick_vsurf_km_per_sec`,
- updates `s%v`, `s%xh(s%i_v,:)`, `s%xh_start(s%i_v,:)`, and `s%v_start`,
- reports max `|v|`, max `|v|/cs`, and a warning when the interior velocity
  exceeds ten times the requested surface kick,
- copies the selected period into `s%rsp_period` or `s%RSP2_period` when
  appropriate.

## Static Background Policy

The current code writes `max_abs_v_div_csound` to the matrix summary but no
longer rejects nonzero background velocities. This matches the RSP workflow,
where the model can be kicked from a nearly static envelope without an arbitrary
velocity/csound cutoff.

Mathematically, the implemented pressure work and turbulent pressure work rows
are static background linearizations. They are the right default for RSP LNA
use. A true nonstatic LNA would need extra first order terms involving the
background velocity, and that is not the target now.

Velocity time centering in the nonlinear Newton residual should not change the
continuous linear eigenproblem assembled here. The LNA does not currently use the
finite dt residual Jacobian, so time centering controls should be allowed as long
as the active continuous equations are otherwise supported.

## Boundary Conditions

Supported or partly supported now:

- Fixed/zero surface velocity paths.
- Surface momentum/pressure boundary variants used by MESA hydro.
- `use_RSP_L_eqn_outer_BC` through the luminosity row:

```text
L(1) = RSP2_Lsurf_factor*4*pi*r(1)^2*c*a*T(1)^4
```

- Surface temperature boundary through `surface_lnT_bc_for_star_LNA`.
- Fixed inner envelope radius/velocity through the current `wrap_r_p1` and
  innermost velocity treatment.

Known boundary condition audit items:

- Verify surface momentum row against `hydro_momentum` for each active outer BC.
- `use_RSP_L_eqn_outer_BC`, or `RSP2_use_L_eqn_at_surface` with `RSP2_flag`,
  selects the surface luminosity boundary shown above.
- Verify the innermost pressure work and momentum face mass for envelope models.
- Full star center boundary conditions remain a later phase.

## Known Open Issues

The items below are not all bugs, but they are places where the code and math
need another pass.

1. The support module remains large.

   The public flow and main row order are isolated in `star_LNA.f90`, and
   turbulence closures are in `star_LNA_turbulence_closures.f90`. A further
   solver/output split is optional and should wait until model validation shows
   that it reduces active maintenance risk.

2. `Eq` heating is intentionally zero in the static gas energy row.

   This is the static background first order perturbation of a quadratic
   velocity gradient heating term. The RSP2 turbulent energy row also omits
   `Eq` from the static first order operator, matching RSP LINA. TDC
   `Eq_div_w` now returns zero for the same reason.

3. TDC MLT turbulent pressure is not present as its own work output component.

   Momentum and static work pressure include it. EOS and RSP2 turbulent pressure
   now share the operator closure with work output, while the TDC face centered
   term still needs a separately named diagnostic.

4. RSP2 cell terms now use a cell pressure scale height.

   The nonlinear RSP2 residual, its pre-solver turbulent velocity estimate,
   eddy viscosity coefficient, and star LNA now form
   `Hp_cell = Peos(k)/(rho(k)*0.5*(grav_face(k) + grav_face(k+1)))`.
   Ordinary damping uses `Hp_cell` and radiative damping uses `Hp_cell^2`.
   Face source and luminosity terms continue to use face scale heights. The
   turbulent pressure row uses
   `etrb = w^2` once through `wrap_etrb_00`; no extra `w^2` was found in that
   row. Numerical comparison against matched RSP/RSP2 model output remains.

5. Mode selection still mirrors RSP's practical growth filter.

   This is intentional for now. If the printed modes are still not matching
   expected Cepheid radial modes, the first response should be row auditing, not
   adding period window controls.

6. Algebraic elimination is documented and checked at source level.

   The Schur complement math and reconstruction minus sign are adjacent to the
   implementation. Numerical conditioning still needs representative model
   output.

7. Existing MESA residuals are not mechanically tied to LNA rows.

   The current approach is explicit continuous equation assembly. That is good,
   but changes in the nonlinear residuals still require a source audit. The
   equation registry, row structure output, and manuscript provide the current
   source code map.

8. Several star branches are intentionally out of scope.

   `check_star_LNA_model` rejects rotation, RTI, mass corrections, user
   `other_*` hooks, velocity drag, and `use_compression_outer_BC`. It also
   rejects RSP2 eddy viscosity with `u_flag` until that nonlinear Riemann term
   is moved into the cell momentum source.

9. The work integral normalization is not final.

   It should be checked against RSP `LINA_work*.data`, including pressure,
   turbulent pressure, and eddy viscosity signs.

## Equation Audit Checklist

Use this checklist before changing selection logic or adding controls.

Density:

- Compare `cell_volume_for_star_LNA` with MESA density/volume closure.
- Check fixed innermost radius behavior for envelope models.
- Confirm the row is algebraic and has no entries in `B`.

Radius:

- Compare with `do1_radius_eqn` in the continuous limit before time centering.
- Confirm `sigma*dlnR = dv/r`.
- Confirm force zero velocity cells map to `dv = 0`.

Momentum:

- Compare `momentum_rhs_for_star_LNA` with `hydro_momentum.get1_momentum_eqn`.
- Confirm sign of gravity and pressure jump.
- Confirm `star_LNA_dm_face` for normal star, RSP, and RSP2 cases.
- Confirm surface momentum/pressure row for each outer BC.
- Confirm artificial viscosity pressure support or intentional omission.
- Confirm TDC MLT turbulent pressure and RSP2 turbulent pressure are not double
  counted.

Energy:

- Compare `energy_rhs_for_star_LNA` with the continuous hydro energy equation.
- Confirm `dL/dm` sign.
- Confirm static pressure work stencil.
- Confirm `dE/dRho*rho` and `Cv*T` partials.
- Confirm TDC/RSP2 turbulent energy inertia.
- Decide whether any first order `Eq` term belongs in a static LNA.

Luminosity:

- Compare surface T and RSP surface luminosity rows with active MESA controls.
- Compare the temperature gradient residual without TDC with the MESA temperature equation.
- Compare RSP2 `Lr + Lc + Lt - L` with `hydro_rsp2`.
- Compare TDC `Lrad + Lconv - L` with `turb:set_TDC_LNA`.

TDC:

- Verify face centered quantities: `mlt_vc`, `Hp`, `gradT`, `L_conv`.
- Confirmed `A = conv_vel/sqrt(2/3)` consistently.
- Confirmed `TDC_alpha_Pt` term placement on B side.
- Confirmed the LNA bypasses the flux limiter derivative branch.
- Confirmed `Eq_div_w` is local, side effect free, and zero in the static
  first order operator.

RSP2:

- Compare `COUPL`, `Lt`, `Ptrb`, `Hp`, and cutoff logic with RSP2 hydro.
- Compare with `rsp_lina` for branch behavior and expected mode ordering.
- Check radiative damping and turbulent pressure work signs.

Solver/output:

- Confirm `sigma = alpha/beta` from `DGGEV`.
- Confirm row/column scaling is undone correctly for eigenvectors.
- Confirm selected modes and raw modes are both visible.
- Confirm kick mode numbers are selected list modes, not raw LAPACK indices.

## Static Source Audit Pass 1

This is the earlier audit record. The current row selection, correlation response,
solver, and viscosity discretization are described above and in the
2026-09-16 source audit below. Historical run results remain as recorded.

This pass compared the current star LNA rows against the nearby MESA hydro source
and RSP LINA without running a model.

### Momentum Source Comparison

MESA source:

- `star/private/hydro_momentum.f90:get1_momentum_eqn`
- `star/private/hydro_momentum.f90:expected_non_HSE_term`
- `star/private/hydro_rsp2.f90:compute_Uq_face`
- `star/private/tdc_hydro.f90:compute_tdc_Uq_face`

MESA's `v_flag` momentum residual is documented in `get1_momentum_eqn` as:

```text
0 = other + grav - RTI_terms - dPtot*area/dm - d_mlt_Pturb*area/dm
```

where:

```text
other = extra_grav - dv/dt + Uq
```

Dropping unsupported hooks (`extra_grav`, RTI, drag) and solving for `dv/dt`
gives:

```text
dv/dt = grav + Uq - (dPtot + d_mlt_Pturb)/(dm/area)
```

This matches the current star LNA `momentum_rhs_for_star_LNA` sign convention.
The current `star_LNA_dm_face` also mirrors MESA's face mass:

```text
k > 1:  dm_face = 0.5*(dm(k) + dm(k-1))
k = 1:  dm_face = 0.5*dm(k)
```

with mass correction logic present in star LNA but currently rejected by
`check_star_LNA_model`.

Open momentum source audit items:

- MESA `get_dPtot_face_info` includes `Pvsc`, RSP2 turbulent pressure,
  `extra_pressure`, and time centering where active. star LNA currently omits
  artificial viscosity and extra pressure by policy. That is RSP, but the
  exact policy should stay explicit.
- MESA's MLT turbulent pressure in momentum uses `mlt_vc_old(k)` in the normal
  helper. star LNA uses the internal TDC `w` when TDC LNA is active. That is the
  desired perturbed TDC behavior, but it differs from the value used by the
  nonlinear residual from the previous timestep. Check the difference with row
  audit output.
- `Uq_face_for_star_LNA` uses `s%dm_bar(k)` inside the Uq helper, matching both
  RSP2 and TDC hydro helpers. The pressure acceleration uses `star_LNA_dm_face`.
  The `Chi` perturbation is reconstructed locally as a static coefficient times
  `delta(v/r)` instead of using cached hydro `Chi_ad`. This split matches the
  static LNA target and should be labeled in the cleanup.

### Radius and Time Centering Comparison

MESA source:

- `star/private/hydro_momentum.f90:do1_radius_eqn`
- `star/private/evolve.f90`, where `using_velocity_time_centering` becomes true
- `star/private/star_LNA.f90:report_star_LNA_setup`

The current star LNA radius row uses the continuous equation:

```text
d lnR/dt = v/r
```

It does not use the finite dt Newton residual or velocity time centering. This is
the right target for a continuous eigenproblem. Time centering controls should not
be a hard blocker for LNA as long as the background quantities have been prepared
and the active continuous equations are supported.

### Energy and Work Source Comparison

MESA source:

- `star/private/hydro_energy.f90:get1_energy_eqn`
- `star/private/hydro_energy.f90:eval_dwork`
- `star/private/hydro_energy.f90:eval1_work`
- `star/private/hydro_energy.f90:eval_simple_PdV_work`

MESA's non-`eps_grav` dE/dt form has the same main sign as star LNA:

```text
esum = -dL_dm + sources + others - d_turbulent_energy_dt - dwork_dm - de_dt
```

Before the P*d(1/rho) time centering branch is active, the MESA `dedt` residual
also contains `-dke_dt - dpe_dt`. Moving these time derivatives to the sigma
side gives the current star LNA structure:

```text
delta[-dL/dm + sources - dwork/dm] = sigma*delta e_eff
```

For RSP2, `e_eff = e + etrb + mechanical`. Since the separate turbulent row is

```text
sigma*delta etrb =
   delta[COUPL - dLt/dm - Ptrb0*dVdt/dm]
```

the total energy row and turbulent energy row can be linearly combined to
recover the gas energy form:

```text
sigma*delta e =
   delta[-dLr/dm - dLc/dm - COUPL - Peos0*dVdt/dm]
```

This is the same structural cancellation used by RSP LINA: `Lt` and turbulent
pressure work live in the total/turbulent split, not as counted twice gas
energy terms.

MESA `eval_dwork` builds:

```text
dwork = [P*A*v]_k - [P*A*v]_{k+1}
```

and its `eval1_work` pressure includes EOS pressure, artificial viscosity, RSP2
turbulent pressure, MLT turbulent pressure, and extra pressure depending on
controls. The current star LNA work row uses the static background limit:

```text
delta dwork = P0*delta(A*v)_k - P0*delta(A*v)_{k+1}
```

That is the correct first order row when `v0 = 0`. It intentionally ignores
`delta(P)*A*v0` terms. If a future LNA supports a genuinely nonstatic
background, this row must be replaced by a full face work perturbation.

Energy source audit items:

- `hydro_energy` includes `Eq_ad` from RSP2 or TDC when those physics are active.
  star LNA currently returns zero from
  `turbulent_viscous_heating_for_star_LNA`. This is correct for the static
  first order operator because eddy viscosity heating is quadratic in
  velocity gradient perturbations, matching the RSP LINA treatment.
- Hydro commit `74842200` also adds total-energy work proportional to `v*Uq`.
  On a static background, radial velocity and `Uq` are first order, so their
  product is second order. The corresponding first-order Star LNA source is
  zero. This is the same order counting that removes `Eq`.
- `hydro_energy` has `others_ad` for diffusion, sedimentation, premixing, phase
  separation, RTI diffusion, drag energy, and mass change heating. star LNA
  omits or rejects those paths. That is acceptable for the current Cepheid
  envelope target but should be part of `check_star_LNA_model`.
- `eval1_work` includes artificial viscosity when `use_Pvsc_art_visc` is true.
  star LNA reports that it ignores artificial viscosity pressure like RSP LINA.
  If exact compatibility with the MESA Newton residual is required, this becomes an
  option or a row term.

### Temperature and Surface Luminosity Source Comparison

MESA source:

- `star/private/hydro_temperature.f90:do1_dlnT_dm_eqn`
- `star/private/hydro_temperature.f90:set_RSP_Lsurf_BC`
- `star/private/hydro_temperature.f90:eval_dlnPdm_qhse`

MESA temperature gradient equation:

```text
resid = delm*dlnPdm_qhse*gradT - lnTdiff
```

with:

```text
lnTdiff = (T(k-1) - T(k))/Tpoint
Tpoint = alfa*T(k) + (1-alfa)*T(k-1)
```

For the logarithmic-pressure form, star LNA uses the same structure through
`temperature_gradient_resid_for_star_LNA` and
`star_LNA_eval_dlnPdm_qhse`. For `dPrad/dm`, static MLT and RSP2 use
`dPrad_dm_resid_for_star_LNA`. Perturbed TDC instead inverts the same transport
relation to obtain the spatial gradient before evaluating the coupled TDC
luminosity and turbulent-velocity relations.

MESA's RSP luminosity outer BC is:

```text
L_1 - RSP2_Lsurf_factor*area*c*a*T_surf^4 = 0
```

scaled by the maximum starting luminosity. star LNA uses the same surface
luminosity condition when `use_rsp_lsurf_row_for_star_LNA` is active, although
the exact sign and scale should remain visible in the row comment after the
cleanup.

Open temperature source audit items:

- MESA `eval_dlnPdm_qhse` time centers EOS pressure when velocity time centering
  is active. star LNA uses the continuous equation without time centering. This
  should not change the eigenproblem, but if the loaded background
  was prepared with time centered quantities, row residual diagnostics can differ.
- MESA `eval_dlnPdm_qhse` uses `mlt_vc_old` from the previous time for MLT turbulent
  pressure in thermodynamic gradients. star LNA uses internal TDC `w` when TDC
  LNA is active. This is physically desirable for perturbed TDC, but it is a
  deliberate difference from the Newton residual.

### TDC Source Comparison

MESA source:

- `turb/public/turb.f90:set_TDC_LNA`
- `star/private/tdc_hydro.f90:compute_tdc_Eq_div_w_face`
- `star/private/tdc_hydro.f90:compute_tdc_Uq_face`

The public wrapper `set_TDC_LNA` is now the central TDC LNA relation. It uses:

```text
A = conv_vel/sqrt(2/3)
L_rad = L0*gradT_actual
L_conv = alpha_C*alpha_eff*alpha_c*rho*T*Cp*4*pi*r^2*A*Y_env
velocity_rhs = S0 - A*DR0 - A^2*D0
velocity_inertia = 2*A + optional alpha_Pt term
```

When `use_dPrad_dm_form_of_T_gradient_eqn` is active,
`actual_gradT_for_star_LNA` first obtains

```text
Lrad = -c*area^2*lambda*dPrad/(dm_bar*kap_floor)
gradT_actual = Lrad/L0
```

from the reconstructed face state and the adjacent cell temperatures. It then
passes this gradient to `set_TDC_LNA`. This preserves both the radiation
transport equation and the single TDC closure row without adding another
unknown.

The alpha_Pt term belongs on the B side because it multiplies
`sigma*delta(1/rho_face)`.

star LNA passes `.false.` for the TDC enthalpy flux limiter.  The limiter is a
nonlinear hydro saturation device; the linear operator should use the
unsaturated `S0` and `L_conv` relation.

Open TDC source audit items:

- `tdc_Eq_div_w_for_star_LNA` is intentionally local and returns zero in the
  static first order operator. The normal TDC hydro source can compute
  `Eq_div_w`, but that term is quadratic in the velocity gradient perturbation,
  so it follows the same omission as RSP LINA's `Eq` term.
- Main at `ee1f16e2` removed the alternate density form TDC eddy viscosity
  control and implementation. Star LNA and nonlinear TDC hydro now both use
  the velocity form. `TDC_alpha_M` controls whether the TDC eddy viscosity
  `Uq`/work term is present, and its length factor is the active `Lambda`.
- The hydro branch harmonic dissipation length is included in the Star LNA TDC
  source, convective luminosity, damping, radiative damping, and eddy-viscosity
  coefficient. GYRE schema 130 still exports the legacy `alpha_MLT` and
  `Hp_face` pair and warns when that differs from MESA's active `Lambda`.
- TDC quantities are face centered, while most base variables are cell centered.
  Every wrapper using `wrap_*_m1`, `wrap_*_00`, or `shift_p1` should be audited
  for face/cell placement before changing mode selection.

### RSP LINA Source Comparison

RSP source:

- `star/private/rsp_lina.f90`
- `star/private/rsp.f90:get_LINA_info`
- `star/private/rsp_build.f90`

RSP uses an ordinary dense eigenproblem:

```text
LLL*X = sigma*X
X = {dR, dU, dT, dw}
```

The star LNA generalized form is needed because it carries algebraic variables
like `lnd`, `L`, and `Hp`.

RSP has the same dense all roots scaling issue at a fixed zone count. In
practice RSP LINA is usually fast because RSP envelope models are much smaller
than normal MESA/star models. The star LNA dense solver is restricted by memory
and cubic eigensolver cost. Large meshes require a targeted, banded, sparse, or
iterative solver.

RSP's mode selector sorts by positive `WI` and applies growth filters. It does
not use Q to classify modes. star LNA intentionally mirrors this behavior. If the
selected star LNA modes are wrong, that points back to row physics or variable
placement before it points to missing Q filters.

RSP work output computes:

```text
QWK   = -pi*dm*Im(conjg(dP)*dV)
QWKPT = -pi*dm*Im(conjg(dPtrb)*dV)
QWKEV = Smolec/Kuhfuss eddy viscosity form
```

star LNA divides each work column by modal kinetic energy. The pressure and
turbulent pressure signs follow RSP. Eddy viscosity and luminosity work have not
been compared with matched RSP output, so their sum is not an independent
growth rate estimate.

### RSP2 Source Audit Details

MESA source:

- `star/private/hydro_rsp2.f90:compute_Y_face`
- `star/private/hydro_rsp2.f90:compute_PII_face`
- `star/private/hydro_rsp2.f90:compute_Source`
- `star/private/hydro_rsp2.f90:compute_D`
- `star/private/hydro_rsp2.f90:compute_Dr`
- `star/private/hydro_rsp2.f90:compute_Lc_terms`
- `star/private/hydro_rsp2.f90:compute_Lt`
- `star/private/star_utils.f90:calc_Ptrb_ad_tw`

RSP2 `Ptrb` in MESA is:

```text
Ptrb = RSP2_alfap*(2/3)*rho*etrb
```

which matches `rsp2_Ptrb_for_star_LNA`. This is also consistent with RSP LINA
if RSP's `ALFAP` is interpreted as the full coefficient multiplying `Et/Vol`.

RSP2 convective source and flux are face/cell mixed. LNA uses the same
hydro AD source factor:

```text
PII_face_LNA = s%PII_ad(k)
Source_k = (w_k + seed)*<PII/Hp>_cell*T_k*(P*QQ/Cp)_k
Lc_face = w_face*area*(x_ALFAC/x_ALFAS)*(T*rho)_face*PII_face
Lt_face = -alpha*alpha_t*(area*rho_face)^2*Hp_face*w_face*detrb/dm_bar
```

STAR LNA computes `Lc`, `Source`, `D`, and `Dr` from this AD state and uses
MESA's `Lr_ad` and `Lt_ad`. The RSP2 `Eq` term is not included in the static
first-order turbulent energy row because it is quadratic in the velocity
gradient perturbation. Face quantities follow
`get_RSP2_alfa_beta_face_weights` and the active RSP2 centering control.
It also zeros `k == 1` and `k == nz`, matching `hydro_rsp2:compute_PII_face`
and the RSP inner boundary convention. It should not apply the
forced nonturbulent cell cutoff internally: the luminosity and turbulent energy
rows decide where `w`, `Lc`, `Lt`, and the source are active, while the hydro
source average can still use the adjacent PII face for the last active turbulent
cell.

Radiative damping comparison:

RSP LINA uses:

```text
D_rad = 4*sigma_SB*(gamma_r/alpha)^2*T^3*(1/rho)^2*Et
        /(Cp*kappa*0.5*(Hp_face(k)^2 + Hp_face(k-1)^2))
```

MESA RSP2 and star LNA now use the native cell pressure and density with the
gravity averaged from the two bounding faces:

```text
grav_cell = 0.5*(G(k)*m_grav(k)/r(k)^2
                 + G(k+1)*m_grav(k+1)/r(k+1)^2)
Hp_cell = Peos(k)/(rho(k)*grav_cell)
Dr = w^2*4*sigma_SB*(RSP2_alfar*x_GAMMAR/alpha)^2
     *T^3/(rho^2*Cp*kappa*Hp_cell^2)
```

This replaces the former square of an averaged face scale height with a directly
collocated cell quantity. The nonlinear residual, automatic differentiation
Jacobian, pre-solver turbulent velocity estimate, eddy viscosity coefficient,
and star LNA operator use the same definition. RSP uses an average of squared
face scale heights as its discrete approximation to the same local factor.

Turbulent pressure comparison:

RSP LINA uses:

```text
PTURB = ALFAP*Et/Vol = ALFAP*rho*Et
```

and differentiates this with respect to radius/volume and turbulent energy. MESA
RSP2 uses the automatic differentiation `etrb = w^2` expression. star LNA also uses
`wrap_etrb_00`, so there is no obvious extra `w^2` in the star LNA turbulent
pressure row itself. If an extra `w^2` remains, the next static place to inspect
is any code that multiplies `Ptrb` by `etrb`, `w`, or `Ptrb_div_etrb` after
`calc_Ptrb_ad_tw`.

### RSP2 Growth Mismatch Follow up

Status: added after the TDC kick growth agreed with the nonlinear KE growth per
cycle while the RSP2 kick growth did not. The first actionable finding is a
growth convention mismatch: `star_LNA` printed the RSP LINA `eta`, while the
nonlinear checks can report exact KE fraction, symmetric `GREKM`, or radius
amplitude growth.

Comparison targets:

```text
logKE_per_cycle = 4*pi*sigma_real/sigma_imag
KE_fractional_expected = exp(logKE_per_cycle) - 1
GREKM_expected = 2*tanh(logKE_per_cycle/2)
amplitude_fractional_expected = exp(logKE_per_cycle/2) - 1
```

For small `logKE_per_cycle`, all three are approximately `logKE_per_cycle`,
`logKE_per_cycle`, and `logKE_per_cycle/2`, respectively. Any nonlinear
comparison must use the same diagnostic
definition, should use the same mode period, and should skip the first few
cycles if the kick excites transient content from other modes.

Static findings from the first follow up pass:

- The RSP2 PII face helper in star LNA must use
  `get_RSP2_alfa_beta_face_weights`, not the generic star face weights. This is
  now the intended implementation.
- The RSP2 PII face helper should zero only the true boundary faces
  `k == 1` and `k == nz`, matching `hydro_rsp2:compute_PII_face`. The forced
  nonturbulent cell cutoff belongs in the `w`, `Lc`, `Lt`, and source callers,
  not inside PII itself, because the last active turbulent cell can use the
  adjacent PII face in the source average.
- The total energy/turbulent energy split still appears structurally consistent:
  the star LNA total energy row carries total luminosity and total pressure
  work, while the separate RSP2 turbulent energy row carries `dLt/dm` and
  turbulent pressure compression work. Combining the two rows should cancel the
  turbulent luminosity and turbulent pressure work pieces in the gas energy
  limit, as in the hydro residual.
- The nonlinear RSP2 utilities already track growth per cycle and smoothed
  growth average diagnostics. Leave that runtime logic unchanged; the LNA
  issue is choosing the matching analytic conversion from `logKE_per_cycle` for
  whatever nonlinear diagnostic is being compared.
- Star LNA period/growth files now keep the existing RSP logarithmic kinetic energy column and
  add exact derived columns for KE fractional growth, `GREKM`, and amplitude
  fractional growth, so nonlinear comparisons no longer depend on interpreting a
  generic `growth` label.
- The current star LNA kick sets only the velocity field. That mirrors the old
  RSP practical kick, but for RSP2 the selected eigenmode also has physical
  `w` and `Hp` components. A velocity only kick can therefore have a larger
  turbulent transient component than TDC, where the extra LNA convection variable
  is not a persistent hydro state.

Next static audit items:

1. Confirm which nonlinear diagnostic is being used in each RSP2 run. In
   `star/rsp2_utils`, the terminal/history label `growth` is radius amplitude
   growth (`delta_R_growth_avg`), not KE growth; compare it to
   `amplitude_fractional_growth_per_period`, not to the KE columns.
2. The implemented `<prefix>_rsp2_term_audit.data` records the active LNA
   closure values for `PII`, `Lc`, `Lt`, source, dissipation, radiative
   damping, turbulent pressure coefficient, pressure work divergence,
   `dLt_dm`, RHS, and inertia. A comparison from a model run still needs the
   corresponding nonlinear hydro values from the same saved model.
3. Add selected mode diagnostics that report the relative size and phase of the
   RSP2 `w` and `Hp` eigenfunction components compared with the velocity
   component. Large `w` or `Hp` components would make a velocity only kick a
   poor clean eigenmode initialization.
4. Audit `Lt_ad` usage. Star LNA still reuses MESA's `s%Lt_ad`; if the growth
   mismatch survives the PII fix, compare the AD partials in `compute_Lt` with
   the LNA row entries and the total energy/turbulent energy cancellation.
5. Audit turbulent pressure work by comparing the LNA static coefficient
   `Ptrb0*d(1/rho)/dt` with `hydro_rsp2:setup_Ptrb_dV_ad` in the zero velocity
   limit.
6. Audit mode selection. RSP2 adds strongly damped turbulent branches, so a
   selected positive frequency root can have the right period range but still be
   a mixed acoustic/turbulent mode. Use raw mode output plus the selected mode
   eigenfunction diagnostics before adding any period window selector.

Checks to run after the static audit:

- Compare `grekm_growth_per_period` to the measured nonlinear `GREKM`/symmetric
  KE growth; compare `ke_fractional_growth_per_period` to
  `(KE_n-KE_{n-1})/KE_{n-1}`; compare `amplitude_fractional_growth_per_period`
  to radius amplitude growth.
- Measure growth after the transient cycles, or compare growth in the first and
  later cycles to identify kick impurity.
- Test velocity only versus a future full eigenfunction RSP2 kick if the
  selected modes have large `w` or `Hp` components.

### TDC Eddy Viscosity Audit Details

MESA source:

- `star/private/tdc_hydro.f90:compute_Chi_cell`
- `star/private/tdc_hydro.f90:compute_Chi_div_w_face`
- `star/private/tdc_hydro.f90:compute_tdc_Eq_div_w_face`
- `star/private/tdc_hydro.f90:compute_tdc_Uq_face`

TDC momentum Uq uses the cell `Chi` difference:

```text
Uq_face = 4*pi*(Chi_{k-1} - Chi_k)/(r_face*dm_bar)
```

The current star LNA `static_eddy_Uq_face_for_star_LNA` no longer calls
`compute_Chi_cell` directly. It reconstructs the static `Chi` perturbation
locally so TDC Uq does not inherit normal hydro/TDC cached convective velocity
partials or nonstatic background velocity gradient terms.

Normal TDC turbulent heating source uses `Chi_div_w_face`:

```text
Eq_div_w = 4*pi*Chi_div_w_face*d_v_div_r/dm_bar
```

The star LNA helper does not reconstruct this term for the matrix. It returns
zero because the source is quadratic in the velocity gradient perturbation.
This avoids generating artificial first order terms from small nonzero
background velocities in an otherwise static LNA.

### Automatic Differentiation Mapping Audit

MESA source:

- `star/private/auto_diff_support.f90:wrap`
- `star/private/auto_diff_support.f90:shift_m1`
- `star/private/auto_diff_support.f90:shift_p1`
- `star/private/star_LNA.f90:ad_index_to_star_LNA_var`

The current mapping from automatic differentiation to LNA variables covers:

```text
i_lnd_m1/00/p1 -> lna_var_lnd at k-1/k/k+1
i_lnT_m1/00/p1 -> lna_var_lnT at k-1/k/k+1
i_w_m1/00/p1   -> lna_var_w   at k-1/k/k+1
i_lnR_m1/00/p1 -> lna_var_lnR at k-1/k/k+1
i_v_m1/00/p1   -> lna_var_v   at k-1/k/k+1
i_L_m1/00/p1   -> lna_var_L   at k-1/k/k+1
i_Hp_m1/00/p1  -> lna_var_Hp  at k-1/k/k+1
```

The current map intentionally does not support rotation, `w_div_wc`, or generic
extra variables. If an active AD expression contains one of those partials, LNA
stops instead of silently dropping it. This is the right failure mode.

TDC internal `w` uses the normal AD `i_w_00` slot. In `tdc_A_for_star_LNA`, the
manual `wrap` call sets:

```text
A = mlt_vc/sqrt(2/3)
dA/dw_00 = 1
```

so TDC and RSP2 both use `lna_var_w`, but only one branch is active in a
given model:

```text
RSP2_flag      -> physical RSP2 w variable
TDC LNA   -> internal LNA A variable
```

No incorrect automatic differentiation index mapping was found in this pass. The implemented row structure
file reports the nonzero count and dominant variable for each row of `A` and
`B`.

### Static Row Audit, 2026-05-05

This pass was source only. No MESA compile or model run was done.

Density:

- MESA residual: `hydro_eqns.f90:do1_density_eqn` solves
  `lnR_actual - log(r_inner^3 + dm/rho/(4*pi/3))/3 = 0`.
- `star_LNA` row: `lnrho + log(4*pi/3*(r_outer^3-r_inner^3)) = const`.
- Status: equivalent first order mass/volume closure. Scaling differs but does
  not change the algebraic row.

Radius:

- MESA residual: `hydro_momentum.f90:do1_radius_eqn` uses
  `dr/r0 = v*dt/r0` with optional time centering in the finite dt solve.
- `star_LNA` row: `delta(v_face/r) = sigma*delta lnR`. For `u_flag`,
  `v_face` is the Riemann contact velocity reconstructed from cell variables.
- Status: correct continuous limit for both velocity grids. The LNA disables
  hydro radius and pressure time centering because it is not linearizing a
  finite Newton step.

Momentum:

- MESA residual: `hydro_momentum.f90:get1_momentum_eqn` builds
  `other + grav - RTI - (dPtot+dPmlt_turb)*A/dm_face`.
- The `u_flag` residual is `hydro_riemann.f90:eval_Riemann_dudt_rhs`, with
  pressure flux, geometry, gravity, and the STAR LNA cell eddy-viscosity source
  divided by `dm(k)`.
- `star_LNA` face row: `grav + Uq - (dPtot+dPmlt_turb)/(dm_face/A)`.
- Status: matches the supported static face and cell velocity branches after
  dropping excluded `other_*`, RTI, mass correction, drag, and rotation terms.
  Surface momentum, fixed pressure, and fixed velocity BCs are represented;
  compression outer BC is explicitly out of scope.

Energy:

- MESA residual: `hydro_energy.f90:get1_energy_eqn` uses
  `-dL/dm + sources + others - dEt/dt - dwork/dm - inertia`.
- `star_LNA` row: `-dL/dm + sources - dwork/dm` on A, with thermal,
  turbulent, and mechanical inertia on B where enabled.
- Status: matches the intended continuous static branch. RTI diffusion, drag
  energy, diffusion/phase/premix "others", mass change heating, and `eps_grav` form
  remain out of scope. Eddy viscosity `Eq` is zero at first order for a static
  background.

Luminosity and surface BCs:

- MESA residuals: `hydro_temperature.f90:do1_dlnT_dm_eqn` and
  `set_RSP_Lsurf_BC`; RSP2 uses `hydro_rsp2.f90:do1_rsp2_L_eqn`.
- `star_LNA` rows: logarithmic-pressure or `dPrad/dm` temperature gradient,
  surface temperature, `use_RSP_L_eqn_outer_BC`, and TDC
  `Lrad+Lconv-L`. RSP2 uses the selected temperature-gradient row with its
  reconstructed `Lr`.
- The `dPrad/dm` row uses

  \[
  \Delta P_{\rm rad,expected}
  =-\frac{\Delta m\,\kappa_fL_{\rm rad}}
  {cA_f^2\lambda_f},
  \qquad
  \Delta P_{\rm rad,actual}
  =\frac{a}{3}(T_{k-1}^4-T_k^4),
  \]

  including face reconstruction, the opacity floor, MLT radiative-luminosity
  split, and the optional radiative flux limiter.
- Status: signs and variable placement match the supported equations. RSP2 L
  rows are unscaled relative to hydro; this is algebraically equivalent.

TDC closure:

- Reference: `turb/public/turb.f90:set_TDC` and `set_TDC_LNA`.
- `star_LNA` row: internal `w`/`A` relation with `A = mlt_vc/sqrt(2/3)`,
  unsaturated convective flux, and the `TDC_alpha_Pt` inertia term on B.
- Status: uses the TDC LNA wrapper and has no normal operation side effects.
  The selected logarithmic-pressure or `dPrad/dm` spatial relation supplies
  `gradT_actual`. The hydro enthalpy flux limiter is intentionally ignored.

RSP2 closure:

- Reference: `hydro_rsp2.f90` luminosity, source, damping, radiative damping,
  turbulent pressure, and turbulent luminosity helpers.
- `star_LNA` rows: selected temperature-gradient closure, turbulent energy
  `w`, and algebraic `Hp`.
- Status: source, damping, and `Lt` signs match the static branch. RSP2 hydro and
  star LNA form the scale height for cell damping terms from cell pressure and
  density with gravity averaged from the bounding faces. The finite dt hydro
  `Eq` source is omitted because it is second order in the static LNA
  perturbation.

Algebraic elimination:

- `build_reduced_star_LNA_problem` solves `Aaa*y = Aad*x` and forms
  `Ared = Add - Ada*y`, `Bred = Bdd - Bda*y`.
- `reconstruct_star_LNA_eigenvectors` restores `x_alg = -y*x_dyn`.
- Status: sign is consistent with `Aaa*x_alg + Aad*x_dyn = 0`.

Work output:

- RSP reference: `rsp_lina.f90:LINA_work*.data` writes
  `-pi*dm*Im(conjg(delta P)*delta V)` divided by modal kinetic energy.
- `star_LNA` now writes pressure, turbulent pressure, eddy viscosity, radiative,
  convective, and turbulent luminosity work divided by modal kinetic energy.
- Status: pressure and turbulent pressure signs match RSP's convention. Luminosity
  and eddy viscosity work are diagnostic extensions and still need model
  comparison against RSP/RSP2 outputs.

### Static Checks After Current Cleanup

Commands run, without compiling or running MESA:

```text
git diff --check -- star/private/star_LNA.f90 star/private/star_LNA_support.f90
git diff --check -- notes/star_LNA_plan.md notes/star_LNA_compression_plan.md
awk line length checks on touched Fortran and markdown
rg non ASCII checks on touched Fortran and markdown
```

All passed in the 2026-05 source only cleanup pass. The 2026-08 cleanup ran
`fortitude check`, repeated `git diff --check` and 132-column checks, and built
the `star` library after explicit user permission. No MESA model was run.

## Current Cleanup Direction

The mechanical cleanup is complete: the top level flow owns one problem object,
the equation map is explicit, row math remains next to assembly, turbulence
closures are separated, and matrix summary output carries equation and closure
audits. Further source reorganization is deferred until MLT, TDC, and
RSP2 model results identify a concrete need.

See `notes/star_LNA_compression_plan.md` for completed and deferred items.

## Full-Model Envelope Domain

Status: implemented in source, pending model validation.

`star_LNA_T_inner` selects an outer, contiguous envelope from a full MESA
model. A nonpositive value keeps the complete model. When a positive value
makes an interior cut, the last active zone `k_inner` satisfies

```text
T(k_inner) <= star_LNA_T_inner
T(k_inner + 1) > star_LNA_T_inner.
```

All zones from the surface through `k_inner` remain in the eigenproblem. If no
cell is hotter than the control, the complete model remains active. Otherwise,
the loaded model below the cut supplies fixed background values, but its
perturbations are zero. The inner conditions are then

```text
delta r(k_inner + 1) = 0
delta L(k_inner + 1) = 0
delta state(k_inner + 1:) = 0.
```

For `v_flag`, `delta v(k_inner + 1) = 0` fixes the inner face velocity. For
`u_flag`, the normal Riemann face reconstruction is retained with perturbations
of excluded core cells set to zero. This is also the state produced when an
envelope eigenfunction is applied as a velocity kick and core cell velocities
are set to zero.

The control is declared in `star_data/private/star_controls_dev.inc`, documented
and defaulted in `star/defaults/controls_dev.defaults`, and transferred through
the controls namelist in `star/private/ctrls_io.f90`.
`setup_star_LNA_var_map` in `star/private/star_LNA_support.f90` selects
`k_inner`. `add_ad_partials_to_matrix` omits columns below the cut while the AD
stencils retain their loaded background values. The equation assembly,
diagnostic audits, eigenfunction output, and work integrals use `map%nz`.
`set_star_LNA_velocity_from_eigenvector` applies an envelope kick, zeroes the
excluded core velocity state, and rebuilds the Riemann face state for `u_flag`.

Implementation checklist:

- [x] Add the default-off `star_LNA_T_inner` control and complete its namelist
  plumbing.
- [x] Set `map%nz` from the first inward crossing of `star_LNA_T_inner`.
- [x] Restrict matrix residual audits and closure audits to `map%nz`.
- [x] Report the selected zone and the temperatures on both sides of the cut.
- [x] Apply velocity kicks in the active envelope and zero excluded core
  velocity state.
- [x] Audit the `star_LNA_T_inner <= 0` path: `map%nz = s%nz`, so the existing
  full-domain matrix and output loops retain their bounds.
- [ ] Verify full-domain periods and growth rates against a pre-control run.
- [ ] Compare periods and growth rates for progressively deeper cuts.

## 2026-08-29 dPrad/TDC correction

The 350-zone Cepheid calculation used perturbed TDC together with
`use_dPrad_dm_form_of_T_gradient_eqn = .true.`. The TDC LNA path nevertheless
called the QHSE branch of `actual_gradT_for_star_LNA`. At the inner envelope
boundary the TDC audit gave a closure radiative luminosity of
`2.79d37 erg/s` for a model luminosity near `1.95d37 erg/s`, producing the
reported normalized luminosity residual of `0.426`.

`actual_gradT_for_star_LNA` now inverts the active `dPrad/dm` equation,
including the opacity floor, face reconstruction, and optional radiative flux limiter,
before calling `set_TDC_LNA`. RSP2 already used
`dPrad_dm_resid_for_star_LNA` with `s% Lr_ad(k)` and required no corresponding
change.

The raw roots from the pre-fix run contain an unstable 11.8309-day mode, but
its reconstructed full-pencil eigenvector failed the `1d-3` residual limit.
The selection diagnostic now reports the first otherwise selectable rejected
root together with its maximum-residual row, zone, and row variable. This
distinguishes a remaining Schur reconstruction error from a physical mode
selection failure.

## 2026-08-29 QHSE face-pressure correction

The subsequent calculation with
`use_dPrad_dm_form_of_T_gradient_eqn = .false.` recovered the unstable
11.8313-day fundamental in the raw roots, but rejected it at the full-pencil
residual limit. It also exposed a separate discretization mismatch in the
ordinary QHSE temperature-gradient path. `star_LNA_eval_dlnPdm_qhse` formed
the face pressure by averaging adjacent cell pressures, while the TDC closure
used reconstructed face pressure in

```math
L_0 = \frac{16\pi ac}{3}\frac{GmT_f^4}{P_f\kappa_f}.
```

Since the spatial relation gives

```math
\nabla_{T,f} = \frac{\Delta T/T_f}
{\Delta m_f[-Gm/(4\pi r_f^4P_f)]},
```

using different pressures makes the radiative luminosity proportional to
`P_interpolated/P_reconstructed`. At zone 298 the two pressures differed by
`2.476d-3`, consistent with the `2.392d-3` maximum normalized luminosity-row
residual. The shared QHSE helper now uses reconstructed face pressure when
`use_face_reconstruction` is active. This matches
`hydro_temperature:eval_dlnPdm_qhse` and keeps the spatial gradient and TDC
radiative coefficient on one face discretization.

## 2026-08-29 algebraic zero-w correction

After the QHSE correction, the background TDC luminosity residual fell to
`9.015d-9`. The unstable 11.8313-day root remained in the raw spectrum, but
its full-pencil residual was `1.557d-2`. The largest error was the energy row
at `k = 348`.

The reduced problem reported 1400 dynamic variables and 700 eliminated
algebraic variables. This classified `w` as dynamic in all 350 zones, although
only zones 78 through 199 use the TDC velocity equation. The other 228 zones
impose

```math
\delta w_k = 0
```

without a time derivative. Leaving these constraints in the generalized
eigenproblem makes its time-derivative matrix singular and retains algebraic
roots in the dense solve.

`partition_star_LNA_indices` now requires a nominally dynamic variable to
have a nonzero row in the full time-derivative matrix. The 228 zero-`w` rows
are therefore eliminated with the density and luminosity constraints. For
this model the next calculation should report 1172 dynamic variables and 928
algebraic variables. This keeps the 122 active TDC `w` equations in the
eigenproblem and removes only the forced zero branch. A fresh MESA install
completed successfully after this change.

The corrected model-200 solve reported 1172 dynamic variables and 928
algebraic variables, as expected for 122 active TDC `w` zones. The unstable
fundamental remained at 11.8313 days with `logKE_per_cycle = 8.209d-2`, but its
full-pencil residual was `2.796d-2`. A later calculation with 118 active TDC
zones reported 1168 dynamic variables and a `4.395d-3` residual for the same
11.83-day mode. The zero-`w` partition is therefore correct, but it is not the
source of the remaining residual.

The next diagnostic must evaluate the eigenpair against the reduced pencil
before reconstructing the algebraic variables. This separates the accuracy of
`DGGEV` from numerical loss in the Schur complement or full eigenvector
reconstruction. Replacing the current AD-generated generalized problem with
RSP LINA's hand-assembled standard operator is not justified by these results.

## 2026-08-29 full-pencil eigenpair refinement

The dense reduced solve remains the spectrum finder. Candidate acoustic roots
with an initial full-pencil residual no larger than `1d-1` are refined against
the scaled full descriptor pencil before applying
`star_LNA_max_eigenvector_residual`.

For `r = (A - sigma*B)*x`, the Newton correction is

```math
(A-\sigma B)\,\delta x-Bx\,\delta\sigma=-r,
\qquad c^\dagger\delta x=0.
```

The largest component of `x` supplies the normalization row. With

```math
y=(A-\sigma B)^{-1}(-r),\qquad
z=(A-\sigma B)^{-1}Bx,
```

the update is

```math
\delta\sigma=-\frac{y_j}{z_j},\qquad
\delta x=y+z\,\delta\sigma.
```

The full Star LNA pencil retains the nearest-neighbor radial stencil, so the
complex matrix is packed for `ZGBTRF` and solved with two right-hand sides by
`ZGBTRS`. The implementation permits four accepted Newton steps and
backtracks a correction that does not lower the scaled residual. The original
eigenpair is retained unless the candidate also lowers the componentwise
residual against the original unscaled equations. The Cepheid development
inlist now requests a `1d-8` final residual.

## Frozen Flux Addition Plan

Status: implemented as a source option, pending MESA model validation.

1. Add `star_LNA_convection_treatment = 'frozen_flux'` as the MESA radial
   radial frozen convective flux control.
2. Keep `tdc_lna_active` restricted to `convection_treatment = 'perturbed'`, so
   `frozen_flux` does not add or solve the internal TDC `w` row.
3. Add a luminosity row for `frozen_flux` in models without RSP2:

```text
Lrad_ad + Lconv0 - L = 0
```

4. Keep the surface row unchanged at `k = 1`, matching the existing star_LNA
   surface treatment.
5. Reject RSP2 with `frozen_flux` until an RSP2 specific frozen definition is
   written down.
6. Update setup reporting to print when `frozen_flux` is selected.
7. Update controls documentation to clarify that "fixed" convective flux
   (`star_LNA_perturb_convective_flux = .false.`) is a local derivative switch,
   while `frozen_flux` is a distinct row choice.
8. Validate after compilation by comparing matrix summaries:
   `frozen_flux` should have no TDC `w` rows, should use luminosity rows in the
   interior, and should show `delta Lconv = 0` in work/diagnostic output.

## 2026-09-16 source and manuscript audit

Source: local `EbF/star_lna` at `f6d606939`. The manuscript and the current
sections above were revised against this source. No MESA build or model run
was performed. The earlier run record remains in the readiness checklist.

### Convection and transport

`star_LNA_support:assemble_luminosity_rows` uses the selected temperature
gradient row for static MLT and RSP2. RSP2 also solves the independent
flux residual `Lr + Lc + Lt - L = 0`. Either
`use_RSP_L_eqn_outer_BC`, or `RSP2_use_L_eqn_at_surface` with `RSP2_flag`,
selects the surface luminosity boundary.

The RSP2 closure reads `s%PII_ad`, `s%gradT_ad`, `s%Lr_ad`, and `s%Lt_ad`.
Its signed face relation is

```math
PII_f=x_S(\Lambda_f/Hp_f)C_{P,f}Y_f.
```

`hydro_rsp2:compute_PII_from_Y` evaluates this relation directly, and
`compute_RSP2_gradT` sets gradT=gradL+Y_face. The flux balance is a separate
algebraic row. TDC calls `set_TDC_LNA` with its enthalpy limiter disabled.
RSP2 damping in
`star_LNA_turbulence_closures:rsp2_damping_for_star_LNA` is proportional to
`w^3`, without the former `w_min^3` subtraction.

With `TDC_use_dynamical_gradL`, both convection closures use

```math
Y=\nabla_T-f\nabla_L,\qquad
f=\frac{(dP/dm)_{\rm actual}}{(dP/dm)_{\rm QHSE}}.
```

`turb_support:get_TDC_dynamical_gradL` supplies the AD factor. It leaves the
neutral gradient unchanged at the surface, without `u_flag` or `v_flag`,
for an undefined ratio, and for a nonpositive ratio at an excised physical
inner boundary. The selected MLT pressure coefficient uses `mlt_vc_old`.
The control defaults to false. Its `get_brunt_B` fallback can multiply the
Ledoux composition contribution by `f` twice; that case remains a restriction.

`star_LNA_support:dPrad_dm_resid_for_star_LNA` uses
`(dm(k-1) + dm(k))/2` at `k = nz` with `R_center > 0`.
`star_LNA_turbulence_closures:actual_gradT_for_star_LNA` still uses
`s%dm_bar(k)` for the TDC transport relation. The manuscript distinguishes
these helpers; agreement at the physical inner boundary remains a check.

### Viscosity and work

For face velocity, the operator uses cell stresses and

```math
\delta U_{q,k}=\frac{4\pi}{r_k\Delta\bar m_k}
 (\delta\mathcal X_{k-1}-\delta\mathcal X_k).
```

For cell velocity, it uses face stresses and

```math
\delta\mathcal X_{f,k}=C_{f,k}
 \left(\frac{\delta u_{k-1}}{r_{{\rm mid},k-1}}
       -\frac{\delta u_k}{r_{{\rm mid},k}}\right),\qquad
\delta U_{q,k}=\frac{4\pi}{r_{{\rm mid},k}\Delta m_k}
 (\delta\mathcal X_{f,k}-\delta\mathcal X_{f,k+1}).
```

The implementations are `static_eddy_Uq_face_for_star_LNA`,
`static_eddy_Chi_face_for_star_LNA`, and `Uq_cell_for_star_LNA` in
`star_LNA_turbulence_closures.f90`. The manuscript gives both coefficients,
their boundary rules, and the optional physical inner-boundary TDC stress.
The cell mixing length includes the next face whenever `k < nz`; at `nz`,
it is `Lambda(nz)/2`.

`star_LNA_support:eddy_viscous_work_for_star_LNA` now evaluates work from the
same acceleration operator. It uses `Re(conjg(du)*dUq)` for cell velocity and
the average of the bounding face contributions for face velocity. Pressure
work uses specific volume, `d(1/rho) = -dlnrho/rho`.

### Solver, domain, and export

The manuscript now includes the full-pencil Newton equations documented in
the 2026-08-29 refinement entry and the `star_LNA_T_inner` domain described
above. The dense reduced solve still computes all eigenpairs. Raw-mode output
is written before refinement; selected modes and kicks use accepted pairs.

Schema-130 column 22, `i_tdc_d_conv_h`, contains

```math
D_{\rm conv}=\frac{L_{\rm conv}H_P}
 {4\pi r^2\rho TC_P(\nabla_T-\nabla_L)}\quad[\mathrm{cm^2\,s^{-1}}].
```

This follows `pulse_gyre:store_tdc_lna_point_env`. The manuscript's former
label `D_conv*Hp` was incorrect. The exporter itself has not changed between
the local and remote STAR LNA tips checked in this audit.

### Remote change not integrated locally

The remote STAR LNA tip was `622075fbf`; `origin/main` was `fd396fd73`.
The local checkout lacks the three September 3--4 main commits. Its remote
branch contains them and has rebased history.

Remote commit `39fd1207a` retains the atmospheric face-to-center temperature
offset when the momentum outer boundary is enabled:

```math
T_{\rm bc}=T_{\rm atm}+\frac{Gm\Delta m}{8\pi r^4}
 \frac{\nabla T}{P_{\rm eos}}.
```

Locally, `surface_lnT_bc_for_star_LNA` sets this offset to zero for
`use_momentum_outer_BC`. The remote commit also adds the selected radiation
pressure floor and its derivatives to `atmosphere_surface_P_bc_for_star_LNA`.
The manuscript labels this correction as pending integration; it is not
presented as local behavior.

### Documentation checks

- [x] Update the manuscript equations and routine references.
- [x] Reconcile the implementation notes and validation summary.
- [x] Rebuild the 15-page PDF without LaTeX warnings and inspect every page.
- [x] Synchronize `notes/star_LNA_manuscript.pdf` and
  `output/pdf/star_LNA_manuscript.pdf`; their SHA-256 checksums match.
- [x] Check the edited TeX and Markdown files for whitespace errors.

## Manuscript equation sequence, 2026-09-18

Reorganized `star_LNA_manuscript.tex` around the rows assembled at local
commit `f6d606939`. Each discrete model equation is followed immediately by
its linearization and the implementing routine. Blue identifies dynamic
rows; rust identifies algebraic rows. Text labels preserve the distinction
in grayscale. Constitutive relations, local closure solves, and diagnostics
remain black because they do not add independent global rows.

```math
\delta F=\sigma\,\delta Q\quad\Longrightarrow\quad
\mathbf A\,\delta\mathbf x=\sigma\mathbf B\,\delta\mathbf x,
\qquad
\delta G=0\quad\Longrightarrow\quad\mathbf B_{G,:}=0.
```

The main sequence is density, radius, momentum, energy, transport, the
additional TDC or RSP2 rows, and boundary or forced-zero replacements.
`star_LNA_support` supplies the row assemblers. The transport alternatives
retain the precedence in `temperature_gradient_resid_for_star_LNA`; the TDC
luminosity closure retains its own global row. The RSP2 scale-height equation
now gives both expected-height branches from
`hydro_rsp2:Hp_face_for_rsp2_eqn`. The assembly and eigensolver follow the
complete row sequence. Longer closure definitions, output conventions,
source references, and the existing validation record follow in appendices.

The source audit remains dated 2026-09-16. This presentation revision does
not integrate the remote surface-boundary correction or change MESA source,
controls, model results, or validation status. No MESA build or model run
was performed.

- [x] Pair each model equation with its linearized matrix row.
- [x] Check the displayed rows and branch choices against their source routines.
- [x] Build the revised PDF without LaTeX warnings.
- [x] Inspect all 23 rendered pages; keep each equation pair, its definitions,
  and its code reference together. Verify 7 dynamic and 10 algebraic pairs.
- [x] Synchronize the notes and output PDF copies and verify their checksums.
- [x] Check the revised TeX and Markdown for whitespace errors and resolve
  all equation and section references.

## Manuscript energy and velocity distinctions, 2026-09-18

The equation layout is retained, with the energy-form distinctions restored
next to the energy row. `energy_eqn_option = 'dedt'` is required by the
current LNA; `check_star_LNA_model` rejects `eps_grav_form_for_energy_eqn`.

For the local-work control, the paired substitutions are

```math
e_{\rm eff}=e+E_t,\qquad
\delta\mathcal W_{\rm local}
=\frac{P_{c,0}}{\Delta m}
 \left(\mathcal A_k\delta v_{f,k}
 -\mathcal A_{k+1}\delta v_{f,k+1}\right).
```

For flux work, they are

```math
e_{\rm eff}=e+E_t+K+\Phi,\qquad
\delta\mathcal W_{\rm flux}
=\frac{P_{f,k,0}\mathcal A_k\delta v_{f,k}
 -P_{f,k+1,0}\mathcal A_{k+1}\delta v_{f,k+1}}{\Delta m}.
```

`dwork_dm_for_star_LNA` and `mechanical_energy_inertia_for_star_LNA`
select these together. The switch is independent of nonlinear velocity
time centering. Static kinetic inertia vanishes, but potential inertia
does not. The manuscript now includes the discrete kinetic and potential
energies from `star_utils:cell_specific_KE_qp` and `cell_specific_PE_qp`,
the turbulent-energy inclusion controls and TDC face-to-cell average,
and the RSP2 gas/total/turbulent energy cancellation.

The uncolored `eps_grav` comparison is explicitly outside the assembled
LNA system. At fixed composition, with EOS pressure and no turbulent terms,

```math
\delta\epsilon_{\rm grav}
=-\sigma\left(\delta e-\frac{P_{\rm eos}}{\rho}\delta\ln\rho\right).
```

This comparison follows `hydro_energy:get1_energy_eqn`,
`hydro_energy:setup_eps_grav`, `eps_grav:do_std_eps_grav`, and the local
controls documentation. The manuscript distinguishes continuous
thermodynamic identities from finite-step discretization and numerical
EOS consistency. It also explains that `eps_grav` is not `-dPhi/dt`.

For `u_flag`, the main text now distinguishes cell momentum, the
reconstructed contact velocity used in radius/work, and the surface
energy pressure from the momentum boundary pressure. The Riemann formulas
and their AD variation follow `hydro_riemann:do1_uface_and_Pface`.
The physical inner boundary is distinguished from an LNA temperature cut.
The face-momentum MLT pressure term is documented as a convection-velocity
coefficient times a cell-density difference.

This pass changes documentation only. It does not add `eps_grav` support,
change source, integrate the remote correction, or validate a model.

- [x] Trace both energy branches and the velocity-dependent coefficients.
- [x] Restore the distinctions in the main text and supporting derivations.
- [x] Rebuild without LaTeX warnings and inspect all 29 rendered pages.
- [x] Synchronize the notes/output PDF copies and check references and whitespace.
  Both copies are identical. All original labels remain, all references resolve,
  and the 17 equation pairs retain 7 dynamic and 10 algebraic variants.

## Upstream synchronization and local-work restoration, 2026-09-18

The live remote and a fresh fetch both gave `origin/EbF/star_lna = 622075fbf`
and `origin/main = fd396fd73`. The local September 1 tip `f6d606939` had
no uncommitted tracked changes. Its development commits were already
represented in the remote's rebased history, including the equivalent
September 1 tip `bacefbf7b`. The large ahead/behind count did not represent
that many independent changes.

The old branch was preserved as
`codex/backup-star-lna-before-upstream-20260918`. Local untracked work was
saved with `git stash push --include-untracked` in stash
`0e3772856f6df49143066595ad8b0d61b8cdb2fb`, labeled
`codex: local work before upstream sync 2026-09-18`. A separate archive of
the local work directories, including their ignored contents, is retained at
`/private/tmp/mesa-star-lna-sync-20260918-ko95cyu7/local-work.tar`.
That directory also contains the before-state manifest and verification record.

With the work saved, `git reset --keep origin/EbF/star_lna` aligned the
branch with the rebased upstream. The stash was applied and retained.
All 3,917 protected entries were verified against the original manifest,
including content hashes, permissions, symlinks, and empty directories.
No conflicts remain. HEAD equals `origin/EbF/star_lna` and contains
`origin/main`. No MESA compilation or model run was performed.

The STAR LNA source difference is exactly the surface-boundary change in
commit `39fd1207a`. It affects `surface_lnT_bc_for_star_LNA` and
`atmosphere_surface_P_bc_for_star_LNA`:

```math
T_{\rm bc}=T_{\rm atm}+\Delta T_0,\qquad
\Delta T_0=\frac{Gm\Delta m}{8\pi r^4}\frac{\nabla T}{P_{\rm eos}}.
```

The temperature offset now remains present with the momentum outer boundary,
and its full AD response enters the LNA surface-temperature row. The optional
momentum pressure floor uses `max(P_bc_raw, a*T_atm^4/3)` and the derivative
of the selected branch. The other STAR LNA implementation files and the
GYRE exporter are unchanged in the comparison.

The manuscript and current validation summary now describe the synchronized
source. Historical audit entries above retain their original dates and status.
The ambiguous work-output statement was clarified: `total_work` is the
arithmetic sum of six diagnostics in `write_one_star_LNA_work`, not yet a
derived and validated discrete energy identity. Pressure work and thermal
luminosity diagnostics need not be independent contributions. This diagnostic
limitation does not itself establish an error in the eigenvalue growth
`4*pi*Re(sigma)/Im(sigma)`. No work-output implementation was changed.

- [x] Preserve old history and local work before changing the branch.
- [x] Synchronize with upstream and restore all protected local entries.
- [x] Verify branch equality, main ancestry, and absence of conflicts.
- [x] Update the boundary equations and distinguish integration from validation.
- [x] Rebuild the 29-page PDF without LaTeX warnings, inspect all changed
  pages, resolve all references, and synchronize the identical notes/output copies.

## 2026-09-19 RSP3 eigenvector residual audit

The 149-zone user case reports residuals 0.94 to 1 for all selected roots.
Printed mode 12 has P=0.7025478 days and logKE/cycle=0.02594198. Treat its
fundamental identification as a hypothesis to check against its mechanical
eigenfunction; do not accept a near-unit residual solely on that basis.

Trace the full componentwise residual, full/reduced scaling, algebraic
elimination and reconstruction. Reproduce the original matrix in an
isolated copy with a local diagnostic `star_LNA` object, keeping production
source uninstrumented. Inspect the worst rows and distinguish physical
operator errors from loss of relative accuracy in very small components.
Any correction must preserve A*x=sigma*B*x and pass a strict residual check;
do not loosen the selection tolerance to hide the error.

Separately confirmed: terminal `mode` starts at zero, but kick controls
start at one. `star_LNA_kick_mode_1=12` selects printed mode 11. Printed
mode 12 requires kick index 13. `star_LNA_mode_for_period=0` independently
chooses the 1.1677-day period. Keep the user's running case unchanged.


Resolved and installed. The isolated case reproduces the near-unit residual
on the same 149-zone model, P=0.70254783858 days. For the fundamental, the
worst relative defects occur in inner Phi rows whose balanced term magnitudes
are about 1d-20 to 1d-19 of the largest row. Its balanced normwise residual is
3.11172d-11. The existing 0.1 componentwise gate prevented any refinement.
This was not evidence that the period or dominant displacement was wrong.

The gate now retains its former componentwise criterion and also permits
refinement when the balanced normwise residual is <=1d-6:

```math
\epsilon_{\rm norm} =
\frac{\max_i |(A_s x_s-\sigma B_s x_s)_i|}
 {\max_i \sum_j (|(A_s)_{ij}|+|\sigma| |(B_s)_{ij}|)|(x_s)_j|}.
```

The additional normwise criterion only decides whether to try refinement.
The final maximum componentwise backward error is still evaluated against
all rows of the original, unscaled matrices. No denominator floor, relaxed
final tolerance, discarded moment rows or changed physical coefficients is
introduced. Existing candidates satisfying the old gate remain eligible.
The banded Newton refinement and its line search are unchanged.

Validation with the installed native solver, isolated copies, eight threads:

| Final tolerance | Accepted modes | Fundamental period (days) | Fundamental residual |
| --- | ---: | ---: | ---: |
| 1d-8 | 15 | 0.7025478385837853 | 2.04160d-9 |
| 1d-10 | 15 | 0.7025478385837856 | 2.34177d-12 |

The fundamental retains logKE/cycle=0.02594197749. Its eigenvalue changes
by about 1.7d-11 relative to the reproduced initial pair. The long-period
candidates preceding it in the pasted table fail the strict final check.
The fundamental is therefore printed mode 0 after this filtering, requiring
kick index 1 and period index 0. Later accepted roots can still include
entropy-moment modes; frequency sorting is not a radial-order classifier.

`star/private/star_LNA_support.f90` contains the refinement gate and the
optional normwise output of `star_LNA_eigenvector_residual_from_vector`.
The control documentation and this manuscript's refinement section were
updated. Fortitude and `git diff --check` pass. The 35-page PDF was rebuilt
and the changed pages inspected. The diagnostic matrix dump was compiled
only into the isolated `pencil` executable, not the production library.
The user's inlist, executable and LNA products were not changed.

Artifacts: `output/review/star_lna_residual_20260919/`, including copied
user outputs, the full matrix, independent Python residual/refinement
analysis, strict native runs, verified numerical results and install logs.
