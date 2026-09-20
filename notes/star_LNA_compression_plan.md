# star_LNA Source Cleanup and Equation Audit

Status: source cleanup record for branch `EbF/star_lna`, based on
`origin/main` commit `4e05d577`. `make star` completed successfully on
2026-08-12. No MESA model has been run with this implementation.

This file retains the row by row design, equations, source comparisons, and
edit sequence used during the cleanup. Completed and deferred items are marked
in the implementation sections. The shorter validation checklist is in
`star_LNA_readiness_checklist.md`.

star LNA linearizes the continuous radial equations used by MESA/star. It does
not reuse the RSP LINA matrix. RSP defines the comparison output, growth
convention, work convention, and velocity kick behavior.

Rejected model options in the cleanup described by this file:

- The first cleanup did not add `u_flag` LNA support. The later hydro-branch
  integration adds the Riemann cell momentum and face reconstruction path.
- Do not add rotation, RTI, or mass correction perturbations.
- Do not add support for user `other_*` hooks.
- Do not add velocity drag perturbations.
- Do not add `use_compression_outer_BC` support.
- The first cleanup did not add the `dPrad/dm` temperature gradient form. The
  later integration adds the full hydro form.

The current restrictions are maintained in `star_LNA_plan.md` and
`merging_star_lna_with_hydro_changes.md`.

## Implemented Layout

The implementation is now split:

- `star/private/star_LNA.f90` keeps the public entry point and main
  equation order.
- `star/private/star_LNA_support.f90` contains the map setup, row
  implementations, closures, dense solver, output, work/kick, and audit helpers.

The solver eliminates algebraic rows and calls LAPACK `DGGEV` for the full
reduced eigensystem. It has no fixed matrix size cutoff. Runtime and storage
scale with the dense reduced matrix, so saved model validation must begin with
a small mesh.

## Source Rationale

Before this cleanup, `star_LNA.f90` contained matrix setup, row assembly,
closures, the eigensolver, and output. The split exposes the equation order and
keeps solver and output code out of the top level row sequence.

The cleanup preserved the implemented equations. Physics changes require a row
audit and saved model comparison.

## Source Requirements

1. Keep the equation list visible.

   `star_LNA.f90` lists the density, radius, momentum, energy, luminosity,
   RSP2, and TDC rows in assembly order.

2. Put math next to code.

   Each row assembler states the equation, linearized form, and terms placed in
   `A` and `B`.

3. Keep helper interfaces direct.

   Helper types remove repeated matrix and map arguments. Row assemblers retain
   the physical terms and their signs.

4. Keep TDC and RSP2 closures explicit.

   Luminosity, pressure, turbulent energy, and work terms remain grouped by
   closure.

5. Keep solver/output separate from physics.

   Dense matrix scaling, algebraic elimination, eigenvalue sorting, output, and
   kicks remain outside the row assemblers.

6. Use row audits before changing mode selection logic.

   A period discrepancy requires identification of the inconsistent row before
   adding a mode filter.

## File Organization

MESA often favors a small number of files over a deep type hierarchy. The
current split keeps the public flow in the main file and the row/solver helpers
in the support file:

```text
module star_lna
  imports
  public entry points and wrappers
    do_star_LNA
    star_LNA_L_conv_ad wrapper

  main equation assembly
    assemble_density_row
    assemble_radius_row
    assemble_momentum_row
    assemble_energy_row
    assemble_luminosity_row
    assemble_rsp2_rows
    assemble_tdc_rows
end module

module star_lna_support
  constants and small types

  setup
    check_star_LNA_model
    setup_star_LNA_var_map
    allocate/free

  shared geometry and AD insertion
    lna_cell_volume
    lna_face_area
    lna_face_mass
    add_ad_to_A
    add_ad_to_B

  physics closures
    pressure closure
    work closure
    luminosity closure
    RSP2 closure
    TDC closure
    eddy viscosity closure

  solver
    scale
    partition
    eliminate algebraic variables
    DGGEV solve
    reconstruct eigenvectors

  mode/output/kick
    raw roots
    selected period/growth
    eigenfunctions
    work files
    initial velocity kick

  diagnostics
    matrix summary
    row residual audit
end module
```

Further splitting, if required after model validation, is restricted to these
module boundaries:

```text
star_LNA_tdc.f90          TDC closure helpers
star_LNA_rsp2.f90         RSP2 closure helpers
star_LNA_solver.f90       dense generalized eigensolver and elimination
star_LNA_output.f90       period/growth, eigenfunction, work, kick
```

No additional split is required by the implemented operator.

## Implemented Top Level Flow

`star_LNA_support.f90` defines the problem object:

```fortran
type star_lna_problem
   type(star_LNA_var_map) :: map
   type(star_LNA_matrix) :: mtx
   integer, allocatable :: eq_id(:)
end type star_lna_problem
```

`do_star_LNA` follows this sequence:

```fortran
subroutine do_star_LNA(s, ierr)
   call set_vars_if_needed(s, 0d0, 'star_LNA', ierr)
   call check_star_LNA_model(s, ierr)
   call setup_star_LNA_problem(s, problem, ierr)
   call report_star_LNA_setup(s, problem%map)
   call assemble_star_LNA_equations(s, problem, ierr)
   call write_star_LNA_matrix_summary(s, problem, ierr)
   call solve_dense_star_LNA(s, problem%map, problem%mtx, ierr)
   call free_star_LNA_problem(problem)
end subroutine do_star_LNA
```

`setup_star_LNA_problem` builds:

- the variable map;
- the equation map; and
- the dense `A` and `B` matrices.

The equation map is important because the row set and the variable set are not
always identical in meaning. For example, `lna_var_L` is an algebraic luminosity
unknown, while the row occupying that slot may be a temperature gradient row, an
RSP surface luminosity row, a TDC luminosity closure, or an RSP2 luminosity
decomposition.

## Equation Registry

The equation identifiers are:

```fortran
integer, parameter :: lna_eq_density = 1
integer, parameter :: lna_eq_radius = 2
integer, parameter :: lna_eq_zero_velocity = 3
integer, parameter :: lna_eq_momentum = 4
integer, parameter :: lna_eq_surface_momentum = 5
integer, parameter :: lna_eq_surface_pressure = 6
integer, parameter :: lna_eq_surface_fixed_velocity = 7
integer, parameter :: lna_eq_energy = 8
integer, parameter :: lna_eq_luminosity_rsp_surface = 9
integer, parameter :: lna_eq_luminosity_rsp2 = 10
integer, parameter :: lna_eq_luminosity_tdc = 11
integer, parameter :: lna_eq_luminosity_frozen_flux = 12
integer, parameter :: lna_eq_surface_temperature = 13
integer, parameter :: lna_eq_temperature_gradient = 14
integer, parameter :: lna_eq_rsp2_turbulent_energy = 15
integer, parameter :: lna_eq_rsp2_zero_w = 16
integer, parameter :: lna_eq_rsp2_Hp = 17
integer, parameter :: lna_eq_tdc_velocity = 18
integer, parameter :: lna_eq_tdc_zero_w = 19
integer, parameter :: lna_eq_mlt_static_temperature_gradient = 20
```

`setup_star_LNA_equation_map` assigns one identifier to every matrix row.
`write_star_LNA_equation_counts` reports the number of rows of each type.

## Row Comment Standard

Row assemblers use this comment format:

```fortran
! Equation:
!   d lnR_k/dt = v_k/r_k
!
! Linearized form:
!   d(v_k/r_k) = sigma*d lnR_k
!
! Matrix:
!   A(row,:) += d(v_k/r_k)
!   B(row,lnR_k) += 1
```

For algebraic rows:

```fortran
! Equation:
!   L_rad_k + L_conv_k - L_k = 0
!
! Linearized form:
!   d(L_rad_k + L_conv_k - L_k) = 0
!
! Matrix:
!   A(row,:) += d(L_rad + L_conv - L)
!   B(row,:) unchanged
```

These comments state the sign and `A` or `B` placement next to the row code.

## Row Cleanup Record

The equations below remain the design reference. Bullets that propose helper
renames or result types record the original cleanup plan; they are not claims
that those names exist in the source. Implementation status is listed under
`Implementation Phases`.

### Density

Keep the current equation:

```text
d ln rho_k + d ln Delta V_k = 0
Delta V_k = (4*pi/3)*(r_k^3 - r_{k+1}^3)
```

Compression changes:

- Rename `cell_volume_for_star_LNA` to `lna_cell_volume_ad`.
- Move it to a geometry section.
- Add a row residual diagnostic:

```text
G = ln(rho) + ln(Delta V) - ln(dm)
```

The constant is not used in the matrix, but it is useful in the audit file.

### Radius

Keep:

```text
d(v_k/r_k) = sigma*d lnR_k
```

Compression changes:

- Make forced zero velocity an explicit row id or row subtype.
- Put the `dv = 0` branch in the radius row, not hidden behind several helpers.

### Momentum

Target equation:

```text
d[grav + Uq - A_face*(Delta Ptot + Delta Pmlt_turb)/dm_face]
   = sigma*d v_face
```

Compression changes:

- Introduce a local `type(lna_momentum_terms)`:

```fortran
type lna_momentum_terms
   type(auto_diff_real_star_order1) :: grav
   type(auto_diff_real_star_order1) :: area
   type(auto_diff_real_star_order1) :: dPtot
   type(auto_diff_real_star_order1) :: dPmlt_turb
   type(auto_diff_real_star_order1) :: Uq
   real(dp) :: dm_face
end type
```

- Build all terms once in `lna_momentum_terms_at_face`.
- Use the same terms in the row residual audit.
- Keep RSP2 `Ptrb` and TDC MLT `P_turb` separated in the term names to avoid
  double counting.
- Add comments showing which term corresponds to `hydro_momentum`.

Audit against:

- `hydro_momentum.get1_momentum_eqn`
- `compute_Uq_face`
- `compute_tdc_Uq_face`
- RSP `rsp_lina` momentum block

### Energy

Target equation:

```text
d[-dL/dm + sources - dwork/dm] = sigma*d e_eff
```

Compression changes:

- Introduce a local `type(lna_energy_terms)`:

```fortran
type lna_energy_terms
   type(auto_diff_real_star_order1) :: dL_dm
   type(auto_diff_real_star_order1) :: sources
   type(auto_diff_real_star_order1) :: dwork_dm
   type(auto_diff_real_star_order1) :: e_eff
   type(auto_diff_real_star_order1) :: mechanical_e_eff
end type
```

- Keep static pressure work explicitly named:

```text
dwork_static_pressure_dm
```

- Do not imply it is the full nonstatic MESA face work linearization.
- Add an audit comment:

```text
For v0 = 0, d(P*A*v)/dm = P0*d(A*v)/dm.
The implemented row assumes that static background limit.
```

- Keep the static `Eq` helper visibly zero and comment why directly in the
  energy source helper. RSP2 should also omit `Eq` from the static first order
  turbulent energy row, matching RSP LINA.

Audit against:

- `hydro_energy`
- RSP `EX`, `EY`, `EU` matrix rows
- Smolec pressure work convention

### Luminosity

Current behavior has too much branching inside one routine. Replace with row ids:

```text
lna_eq_luminosity_surface_T
lna_eq_luminosity_rsp_surface
lna_eq_luminosity_rsp2
lna_eq_luminosity_tdc
lna_eq_luminosity_tgrad
```

Each row should call one closure and insert one AD residual:

```fortran
select case (eq_id)
case (lna_eq_luminosity_rsp2)
   call lna_rsp2_luminosity_closure(s, k, lum)
   resid = lum%Lr + lum%Lc + lum%Lt - wrap_L_00(s,k)
case (lna_eq_luminosity_tdc)
   call lna_tdc_closure(s, k, tdc)
   resid = tdc%L_rad + tdc%L_conv - wrap_L_00(s,k)
...
end select
```

The RSP surface luminosity row should be labeled:

```text
L_1 - RSP2_Lsurf_factor*4*pi*r_1^2*c*a*T_1^4 = 0
```

or the sign used in code, but the equation and code must match visibly.

### RSP2

Group RSP2 terms in one closure:

```fortran
type lna_rsp2_terms
   type(auto_diff_real_star_order1) :: Lr
   type(auto_diff_real_star_order1) :: Lc
   type(auto_diff_real_star_order1) :: Lt
   type(auto_diff_real_star_order1) :: COUPL
   type(auto_diff_real_star_order1) :: Ptrb
   type(auto_diff_real_star_order1) :: Hp_expected
   logical :: force_zero_w
end type
```

Rows:

```text
Lr + Lc + Lt - L = 0
d[COUPL - dLt/dm - Ptrb0*dVdt/dm] = sigma*d etrb
Hp_expected - Hp = 0
```

Compression changes:

- Build and reuse the same `Ptrb` expression for momentum, work output, and
  turbulent pressure work.
- Keep radiative damping terms next to the RSP2 turbulent source closure.
- Add explicit comments for the `RSP2_nz_div_IBOTOM` cutoff.

Audit against:

- `hydro_rsp2.f90`
- `rsp_lina.f90`
- Smolec/Moskalik equations for turbulent pressure, turbulent flux, and
  radiative damping.

### TDC

Keep the public wrapper in `turb` but reduce the MESA/star wrappers.

Define one MESA/star result type:

```fortran
type lna_tdc_terms
   type(auto_diff_real_star_order1) :: luminosity_resid
   type(auto_diff_real_star_order1) :: velocity_rhs
   type(auto_diff_real_star_order1) :: velocity_inertia
   type(auto_diff_real_star_order1) :: L_rad
   type(auto_diff_real_star_order1) :: L_conv
   type(auto_diff_real_star_order1) :: A
   type(auto_diff_real_star_order1) :: Eq_div_w
   logical :: force_zero_w
end type
```

Then replace:

```text
tdc_luminosity_resid_for_star_LNA
tdc_luminosity_terms_for_star_LNA
tdc_relation_for_star_LNA
```

with:

```fortran
call lna_tdc_terms_at_face(s, k, tdc, ierr)
```

Rows:

```text
L_rad + L_conv - L = 0
d velocity_rhs = sigma*d velocity_inertia
```

Required math comment above `lna_tdc_terms_at_face`:

```text
A = conv_vel/sqrt(2/3)
Y = gradT_actual - gradL
Y_env = Y*Gamma/(1+Gamma) when enabled and Y > 0
L_rad = L0*gradT_actual
L_conv = alpha_C*alpha*alpha_c*rho*T*Cp*4*pi*r^2*A*Y_env
velocity_rhs = S0 - A*DR0 - A^2*D0
velocity_inertia = 2*A + optional alpha_Pt compression term
```

Audit against:

- `turb/public/turb.f90:set_TDC_LNA`
- `turb/private/tdc_support.f90`
- `star/private/tdc_hydro.f90`
- `star/private/turb_info.f90`

Important: the TDC row must remain internal to the LNA. It should not change
normal TDC evolution or add a new regular star variable.

## Shared AD Insertion

Current helpers:

```fortran
add_ad_partials_to_A
add_ad_partials_to_B
ad_index_to_star_LNA_var
```

Keep one path for AD insertion, but make the mapping clearer:

```text
AD index -> zone offset -> LNA var id -> matrix column
```

Add a debug mode that can print, for one requested row:

```text
row id
equation id
nonzero A columns with variable names and zone offsets
nonzero B columns with variable names and zone offsets
```

This is more useful than a global dense matrix dump when checking row mistakes.

## Solver Section

Keep the existing algorithm:

1. Row scale and column scale.
2. Partition dynamic/algebraic variables.
3. Require algebraic rows have zero B entries.
4. Eliminate algebraic variables.
5. Solve reduced `Ared*x = sigma*Bred*x` with `DGGEV`.
6. Reconstruct algebraic components with the minus sign.
7. Undo eigenvector column scaling.

Rename variables for clarity:

```text
Aaa -> A_alg_alg
A_ad -> A_alg_dyn
A_da -> A_dyn_alg
B_da -> B_dyn_alg
alg_from_dyn -> inv_A_alg_alg_times_A_alg_dyn
```

Put this math above the reduction code:

```text
A_alg_alg*x_alg + A_alg_dyn*x_dyn = 0
x_alg = -inv(A_alg_alg)*A_alg_dyn*x_dyn
```

Then the code line:

```fortran
vr_full(alg_idx,j) = -matmul(inv_A_alg_alg_times_A_alg_dyn, vr_red(:,j))
```

will be self checking.

## Output Section

Group outputs into four routines:

```fortran
call lna_write_raw_roots
call lna_select_modes_rsp_style
call lna_write_period_growth
call lna_write_mode_outputs
```

Do not hide raw roots. Keep:

```text
<prefix>_raw_modes.data
```

because it is the best way to distinguish slow roots from the selected list.

Keep the RSP convention:

```text
growth = 4*pi*sigma_re/sigma_im
period_days = 2*pi/sigma_im/86400
```

Do not add period window controls as the primary fix for wrong periods. If the
selected list is wrong, audit the rows first.

## Work Integral Plan

Unify the work output pressure terms with the matrix physics.

Define:

```fortran
type lna_pressure_terms
   type(auto_diff_real_star_order1) :: Peos
   type(auto_diff_real_star_order1) :: Ptrb_rsp2
   type(auto_diff_real_star_order1) :: Ptrb_mlt_tdc
   type(auto_diff_real_star_order1) :: P_total_for_momentum
   type(auto_diff_real_star_order1) :: P_static_for_work
end type
```

Then the work writer can output:

```text
gas_pressure_work
rsp2_turb_pressure_work
tdc_mlt_turb_pressure_work
eddy_visc_work
rad_lum_work
conv_lum_work
turb_lum_work
```

This will make it obvious whether TDC turbulent pressure is being included in the
same way in momentum and diagnostics.

Normalize and sign check against RSP:

```text
RSP: QWK   = -pi*dm*Im(conjg(dP)*dV)
RSP: QWKPT = -pi*dm*Im(conjg(dPtrb)*dV)
RSP: QWKEV = eddy viscosity formula from Smolec Appendix C
```

## Row Audit Diagnostics

Add a developer diagnostic file:

```text
<prefix>_row_audit.data
```

For each row type, write:

```text
row_type
zone
background_residual
row_scale
num_A_nonzero
num_B_nonzero
max_abs_A
max_abs_B
dominant_A_variable
dominant_B_variable
```

For selected zones, also write nonzero entries:

```text
zone_offset
var_name
A_value
B_value
```

This should be controlled by a development logical, or piggyback initially on
`star_LNA_write_matrix_summary` until the control list is cleaned up.

## Equation Audit Against MESA and RSP

Audit each equation before changing physics:

### Density

- Code: `lna_cell_volume_ad`.
- MESA comparison: density/volume closure.
- RSP comparison: `DVR`, `DVRM`, `Vol`.
- Checks: sign, inner boundary, algebraic row has no B entries.

### Radius

- Code: radius row.
- MESA comparison: `do1_radius_eqn` in continuous limit.
- RSP comparison: velocity definition row.
- Checks: `sigma*dlnR = dv/r`.

### Momentum

- Code: momentum terms closure.
- MESA comparison: `hydro_momentum`.
- RSP comparison: momentum matrix block `MX`, `MY`, `MU`, `MZ`.
- Checks:
  - pressure jump sign,
  - gravitational derivative sign,
  - face mass,
  - surface boundary,
  - eddy viscosity acceleration,
  - RSP2/MLT turbulent pressure separation.

### Energy

- Code: energy terms closure.
- MESA comparison: `hydro_energy`.
- RSP comparison: energy block `EX`, `EY`, `EU`, `EZ`.
- Checks:
  - `dL/dm` sign,
  - static pressure work sign,
  - `Cv*T` and `dE/dRho*rho`,
  - MESA `dedt` kinetic/potential inertia when active,
  - turbulent inertia,
  - source terms,
  - `Eq` is omitted from the static first order operator.

### Luminosity

- Code: luminosity row dispatch.
- MESA comparison: `hydro_temperature`, TDC/MLT temperature gradient setup.
- RSP comparison: radiative luminosity derivatives and surface luminosity row.
- Checks:
  - surface T BC,
  - `use_RSP_L_eqn_outer_BC`,
  - TDC `Lrad + Lconv`,
  - RSP2 `Lr + Lc + Lt`,
  - face/cell placement.

### TDC

- Code: `set_TDC_LNA` and MESA/star TDC closure.
- MESA comparison: `tdc_hydro`, `turb_info`, `tdc_support`.
- RSP comparison: compare the Kuhfuss source and damping structure directly.
- Checks:
  - `A` vs `conv_vel` conversion,
  - face centered `rho`, `T`, `P`, `Hp`, `gradT`,
  - `TDC_alpha_Pt` term on B side,
  - enthalpy flux limiter is ignored in LNA,
  - `Eq_div_w` source is explicitly zero in the static first order LNA,
  - no normal operation side effects.

### RSP2

- Code: RSP2 closure.
- MESA comparison: `hydro_rsp2`.
- RSP comparison: `rsp_lina`.
- Checks:
  - `Ptrb = alpha_p*(2/3)*rho*etrb`,
  - radiative damping factor and averaging,
  - turbulent flux `Lt`,
  - convective flux `Lc` with unlimited `PII_face_LNA`,
  - `COUPL` reconstructed locally with unlimited source terms,
  - `Eq` omitted from the static first order linearization,
  - forced nonturbulent cells,
  - `Hp` row.

## Implementation Phases

### Phase 0: Documentation and Comments

- [x] Preserve the implemented equations.
- [x] Add row math comments above existing row assemblers.
- [x] Add this equation audit and the implementation notes.
- [x] Compile only after explicit user permission.

Phase 0 progress:

- Added equation comments above density, radius, momentum, energy, luminosity,
  RSP2 `w`, RSP2 `Hp`, and TDC internal-`w` row assembly in
  `star/private/star_LNA_support.f90`.
- Added the Schur complement algebraic elimination math above
  `build_reduced_star_LNA_problem`.
- Added `<prefix>_row_structure.data` behind the existing matrix summary output.
  It lists each full matrix row, A/B nonzero counts, and the dominant A/B
  variable/zone so row coupling can be inspected without dumping the full dense
  matrix.
- Expanded `notes/star_LNA_plan.md` with a static source audit comparing the
  current rows against `hydro_momentum`, `hydro_energy`,
  `hydro_temperature`, `tdc_hydro`, and `rsp_lina`.
- TDC passes `.false.` for its enthalpy flux limiter to `set_TDC_LNA`.
- RSP2 reuses the hydro AD correlation and source in its linear operator.
- Matched the RSP2 `PII_face` boundary behavior by forcing the local LNA helper
  to zero both `k == 1` and `k == nz`, as `hydro_rsp2:compute_PII_face` does.
- Confirmed the RSP2 turbulent energy row omits `Eq` in the static first order
  operator, matching RSP LINA and the quadratic velocity gradient argument.
- Changed eddy viscosity acceleration assembly to reconstruct static `Chi`
  locally for RSP2 and TDC instead of reusing cached hydro `Chi_ad` or normal
  TDC `compute_Chi_cell`.

### Phase 1: Problem Setup Cleanup

- [x] Introduce a `star_lna_problem` type.
- [x] Move variable map and equation map setup into one section.
- [x] Keep the existing matrix indices and equations.
- [x] Add equation identifiers, row names, and equation count output.

### Phase 2: Row Closure Extraction

- Create `lna_momentum_terms`, `lna_energy_terms`, `lna_rsp2_terms`, and
  `lna_tdc_terms`.
- Replace thin wrappers with one closure per physical row family.
- Keep all equations and signs unchanged.
- Add row audit residual calculations using the same closures.

Status: partially implemented. Pressure components are shared by the operator
and work output. RSP2-term and TDC face audit files call the active closure
routines. The proposed result types were not added because they did not remove
enough repeated row logic.

### Phase 3: Solver Naming Cleanup

- Rename algebraic elimination block variables.
- Put the Schur complement equation above the code.
- Keep `DGGEV` and scaling unchanged.
- Add checks that B is zero on algebraic rows with row/equation names in the
  error output.

Status: the Schur complement equations and reconstruction sign are documented
next to the implementation. Algebraic rows are checked for nonzero `B` entries.
Broader local variable renaming is deferred.

### Phase 4: Output and Work Cleanup

- Group raw roots, selected modes, eigenfunctions, work, and kick routines.
- Unify pressure term closures between rows and work output.
- Split TDC MLT turbulent pressure from RSP2 turbulent pressure in diagnostics.
- Keep RSP frequency sorting and eta conventions, with explicit frequency
  and eta selection controls.

Status: output routines remain grouped in `star_LNA_support.f90`. EOS and RSP2
pressure components are shared with work output. Growth output is labeled
`logKE_per_cycle`; legacy `*_eta` control names remain for inlist compatibility.
Eddy viscosity and luminosity work still require comparison with matched RSP
output.

### Phase 5: Saved Model Physics Validation

The following comparisons have not been run:

- Check momentum pressure signs and face masses against MESA hydro.
- Check energy work and `dL/dm` signs against MESA and RSP.
- Check TDC face centered placement.
- Check RSP2 radiative damping and turbulent pressure averaging.
- Check `TDC_alpha_Pt` term by perturbing only density in the TDC wrapper.
- Check work integral normalization against RSP `LINA_work*.data`.

## Cleanup Results

Source results:

- [x] `star/private/star_LNA.f90` has a main equation assembly
  section.
- [x] Every row assembler has math comments above it.
- [ ] TDC and RSP2 closures are returned through named result types. These types
  are deferred because the existing closure routines already provide the audit
  values without a second data structure.
- [x] The same EOS and RSP2 pressure closures feed matrix rows and work
  diagnostics.
- [x] Matrix summary reports row and equation counts.
- [x] A row structure file identifies dominant `A` and `B` entries by variable.
- [x] Raw positive frequency roots are written.
- [x] Selected period/growth output uses the RSP convention and documents
  its frequency/eta selection filters.
- [x] Kick mode numbers start at one and support three components.
- [x] Unsupported branches fail in `check_star_LNA_model`: rotation, RTI, mass
  corrections, user `other_*` hooks, velocity drag, and
  `use_compression_outer_BC`. Later work supports `u_flag` and `dPrad/dm`, with
  the combinations listed in `star_LNA_plan.md`.
- [x] Hydro `v_flag` is not required for analysis; star LNA introduces a face velocity
  perturbation variable for static backgrounds without hydro.
- [x] `make star` completed on 2026-08-12 after explicit user permission.
- [ ] MLT, TDC, and RSP2 saved model comparisons.

## First Edit Set

The first edit set was mechanical:

1. Add equation ids and equation names.
2. Add `star_lna_problem`.
3. Move map allocation/free into problem allocation/free.
4. Add row math comments without changing formulas.
5. Add row/equation counts to the matrix summary.

This edit set changed source organization and audit output. It did not tune mode
periods or selection thresholds.

Status after the 2026-05-05 cleanup pass:

- Added explicit module sections for setup, matrix assembly, AD mapping, dense
  solve/elimination, mode output, work/kick, row audit diagnostics, and utilities.
- Split dense LAPACK solve from solution handling in
  `finish_dense_star_LNA_solution`.
- Split the implementation into `star_LNA.f90` for the public flow and
  main equation order, with `star_LNA_support.f90` for helper
  implementation.
- Added the missing RSP work normalization by modal kinetic energy.
- Added the static row audit to `notes/star_LNA_plan.md`.

All items in this edit set are implemented.

## Deferred Momentum Edit Set

After the first patch is stable:

1. Add `lna_momentum_terms`.
2. Replace separate pressure/area/gravity/Uq helpers with one momentum closure.
3. Add row audit output for momentum rows.
4. Compare the closure names directly against `hydro_momentum`.

This is the likely place to find sign, face mass, or pressure placement mistakes
that could shift acoustic periods.

## Deferred Energy Edit Set

Then:

1. Add `lna_energy_terms`.
2. Make static pressure work explicitly named.
3. Add the energy row audit.
4. Verify `dL/dm`, `dwork/dm`, and `e_eff` signs.

This is the likely place to find growth rate errors.

## Deferred TDC Edit Set

Then:

1. Add `lna_tdc_terms`.
2. Reduce MESA/star TDC wrappers to one closure call.
3. Add a TDC row audit for `velocity_rhs` and `velocity_inertia`.
4. Check face centered placement against TDC/MLT source routines.
5. Keep the LNA TDC closure on the unsaturated relation regardless of
   `use_TDC_enthalpy_flux_limiter`.
6. Preserve the distinction between the old broad `frozen` diagnostic and the
   explicit `frozen_flux` row. `frozen_flux` should remain a no-`w` TDC
   comparison path with `Lrad_ad + Lconv0 - L = 0`. `perturbed` owns the
   internal TDC `w` row and the full local TDC closure. Keep `frozen_flux`
   unavailable for RSP2 until its frozen flux equation is defined.

This is the likely place to find remaining TDC specific period or growth issues.

## Deferred RSP2 Edit Set

Then:

1. Add `lna_rsp2_terms`.
2. Unify RSP2 pressure and luminosity terms across rows and work output.
3. Recheck radiative damping and turbulent pressure work.
4. Compare a source row map to `rsp_lina`.

This is the likely place to find RSP2 specific growth rate errors.
