# Current Branch Changes

Branch: `EbF/improve_split_merge_amr_mesh`

Last updated: 2026-08-26

This is the authoritative catalog of active changes on this branch relative to
`origin/main`. Update it when implementation status, controls, or PPISN test
settings change. Topic notes retain the detailed derivations and audits.

## 1. Split/Merge AMR Metric Zoning

Status: committed metric implementation, including the shared `MaxLong`
split and prospective-merge limit.

The branch adds opt-in metric zoning for split/merge AMR with cell-centered
`u_flag` hydrodynamics. It combines normalized local `logR` and `logtau`
increments:

```text
dxi_k = w_R*abs(dlnR_k)/Delta_lnR
      + w_tau*abs(dlntau_k)/Delta_lntau,

dxi_target = (w_R + w_tau)/N_baseline.
```

The existing split and merge ratios are retained:

```text
split ratio = dxi_k/dxi_target,
merge ratio = max(dxi_target/dxi_k, dq_min/dq_k).
```

The metric zoning path rejects a merge when the summed metric of the proposed
pair would exceed the existing split threshold:

```text
dxi_pair = dxi_i + dxi_ip,
reject merge if dxi_pair/dxi_target > split_merge_amr_MaxLong.
```

Using `split_merge_amr_MaxLong` for both decisions prevents a permitted merge
from producing a cell that must immediately split again. A trace mode reports
metric ranges, components, selected cells, and guard rejections.

The mode is disabled by default and ignored for `v_flag`. Legacy split/merge
zoning remains unchanged when it is disabled or when both weights are
nonpositive.

Controls:

1. `split_merge_amr_use_metric_zoning_for_u_flag = .false.`
2. `split_merge_amr_metric_logR_weight = 1d0`
3. `split_merge_amr_metric_logtau_weight = 1d0`
4. `split_merge_amr_metric_min_delta_lnR = 1d-9`
5. `split_merge_amr_metric_min_delta_lntau = 1d-9`

Primary files:

- `star/private/adjust_mesh_split_merge.f90`
- `star/defaults/controls.defaults`
- `star/defaults/controls_dev.defaults`
- `star_data/private/star_controls_dev.inc`
- `star/private/ctrls_io.f90`

## 2. Conservative `u_flag` Split/Merge Remap

Status: committed through checkpoint `a7e53425c`.

The `u_flag` merge now preserves mass-weighted momentum:

```text
u_new = (dm_1*u_1 + dm_2*u_2)/(dm_1 + dm_2).
```

The resolved kinetic energy removed by that projection is

```text
Delta_KE = 0.5*dm_1*dm_2*(u_1 - u_2)^2/(dm_1 + dm_2),
```

and is added to the merged internal energy. The merge therefore preserves
momentum and the sum of internal plus resolved kinetic energy.

For a split, a limited linear slope constructs two child velocities. Their
mass-weighted mean equals the parent velocity. The added resolved kinetic
energy is removed from child internal energy. Let `e_R*` and `e_L*` be the
provisional child internal energies and let

```text
e_floor = min(e_iR, e_iC, e_iL)
```

be the minimum pre-split internal energy in the three-cell reconstruction
stencil. The allowed specific kinetic-energy correction is

```text
Delta_e_max = max(0, min(e_R*, e_L*) - e_floor).
```

If the requested correction `Delta_e_KE` exceeds this margin, the child
velocity contrast is scaled by

```text
theta = sqrt(Delta_e_max/Delta_e_KE).
```

The resulting correction preserves momentum and internal plus resolved
kinetic energy without introducing a new internal-energy minimum. When the
margin is zero, both children retain the parent velocity.

This limiter was added after the remesh following accepted PPISN model 4617
created adjacent 245 K and 283 K children. Their `logT` values, 2.389 and
2.452, were below the Ferguson opacity-table minimum `logT = 2.7`, so
`kap_get` returned `ierr = -1` before model 4618 entered the solver. Profile
4600 shows candidate splits for which the velocity reconstruction would have
removed 95--99.6 percent of the parent internal energy. The limiter prevents
that remap-created thermodynamic minimum.

The child masses are closed explicitly with `dm_R = dm_parent - dm_L` before
the momentum and energy projection.

Primary file:

- `star/private/adjust_mesh_split_merge.f90`

## 3. Weak-Gravity MLT/TDC Mixing Length

Status: committed.

The branch adds one regular control:

```text
harmonic_dissipation_length_beta = 0d0
```

Positive values enable the harmonic mixing length

```text
Lambda_0 = alpha_mlt*Hp,
1/Lambda = 1/Lambda_0 + 1/(beta*r).
```

This has the limits

```text
Lambda -> alpha_mlt*Hp  when alpha_mlt*Hp << beta*r,
Lambda -> beta*r        when alpha_mlt*Hp >> beta*r.
```

When enabled, MLT and TDC construct `Lambda` from the raw HSE pressure scale
height `Hp = P/(rho*g)` and do not apply `alt_scale_height_flag` to that mixing
length. The thermodynamic pressure scale height itself is not replaced.

The same `Lambda` is used consistently by ordinary MLT, TDC source and damping
terms, convective flux, `D = v_conv*Lambda/3`, stored mixing length, derived
convective velocity, small-convection-region redo logic, and TDC `Chi` and
matching `Eq`. Thermohaline mixing, semiconvection, RSP, and RSP2 retain their
legacy scales.

Primary files:

- `star/private/star_utils.f90`
- `star/private/reconstructed_face_support.f90`
- `star/private/turb_support.f90`
- `star/private/turb_info.f90`
- `star/private/mix_info.f90`
- `star/private/tdc_hydro.f90`
- `turb/public/turb.f90`
- `turb/private/mlt.f90`
- `turb/private/tdc.f90`
- `turb/private/tdc_support.f90`
- `turb/test/src/test_turb.f90`

## 4. Configurable EOS Density Floor

Status: committed, including the explicit PPISN inlist value.

The previous hard-coded rejection below `logRho = -25` is replaced by

```text
min_logRho_for_eos = -30d0.
```

EOS calls below the selected floor return an error and request retry. Values
below `-30` are rejected when controls are read because sufficiently low
density can overflow high-order density derivatives. The PPISN worktree also
sets the control explicitly to `-30d0`.

Primary files:

- `star/defaults/controls.defaults`
- `star_data/private/star_controls.inc`
- `star/private/ctrls_io.f90`
- `star/private/eos_support.f90`
- `star/test_suite/ppisn/inlist_ppisn`

## 5. TDC `u_flag` Surface Traction

Status: committed in checkpoint `16433fecd`.

The outer TDC eddy-viscosity face is explicitly stress free for `u_flag`:

```text
Chi_1 = 0,
Uq_1 proportional to Chi_1 - Chi_2 = -Chi_2.
```

Previously, the missing exterior cell was represented by zero-valued `u` and
`r` wrappers, after which the zero radius was replaced by `1 cm`. This produced
the artificial surface strain `-u_1/r_1`. Zeroing only `Chi_1` removes the
exterior traction while retaining the force from the surface cell's inner
turbulent face. The `v_flag` momentum path already supplies an explicit zero
exterior stress and is unchanged.

Primary files:

- `star/private/tdc_hydro.f90`
- `docs/source/changelog.rst`

## 6. PPISN Configuration

Status: synchronized at branch checkpoint `b0d4e8bff`; static checks are
complete, but the updated models have not been compiled or run.

The maintained pulse implementation is shared by the test-suite case and the
two development cases in `star/dev_cases_test_TDC/ppisn_models`. Its current
configuration is:

1. Metric split/merge AMR uses `logR` weight `1`, `logtau` weight `0.5`,
   `split_merge_amr_nz_baseline = 2000`, and
   `split_merge_amr_MaxLong = 1.25d0`. The same `MaxLong` limit rejects a
   prospective metric merge that would immediately require splitting; the
   separate merge-guard control has been removed.
2. Direct ejecta removal is the only enabled removal algorithm. A contiguous
   surface layer must lie beyond `1d4 Rsun`, have
   `u > 4*vesc`, and have positive total specific energy. `k_keep` is the first
   retained cell, so ejecta sums and cuts consistently use `1:k_keep-1`.
3. Direct removal stores `T(k_keep)`, the updated `tau_factor`, and
   `max[u(k_keep),2*vesc(R_limit)]` as restart state. The velocity is capped by
   `x_ctrl(17) = 5d4 km/s`. The temporary boundary is cleared before the final
   `star_relax_to_star_cut` reconstruction.
4. Riemann-hydro timesteps are limited to `1d4 s` both before and during an
   active pulse. After the first Riemann phase, the hydro-off `v_flag` path is
   capped at `0.5 yr` directly in `run_star_extras`; the obsolete
   `inlist_after_first_pulse` has been removed.
5. `restore_mesh_on_retry = .false.` in both phase inlists. The separate
   `split_merge_amr_avoid_repeated_remesh = .true.` setting prevents one cell
   from being remeshed repeatedly within a single AMR pass; it does not restore
   the pre-remesh grid after a retry.
6. TDC uses `mixing_length_alpha = 2`, `include_mlt_corr_to_TDC = .false.`,
   `TDC_alpha_M = 0`, explicit `mlt_vc` in the momentum equation, turbulent
   energy in the energy equation, and `harmonic_dissipation_length_beta = 1`.
   The shock-convection threshold is
   `max_abs_du_div_cs_for_convection = 0.03d0`.
7. Superadiabatic reduction and its turnover limiter are selected globally in
   `inlist_ppisn`; `superad_reduction_max_logT = 7d0` disables the reduction
   at and above a start-of-step temperature of `1d7 K`.
8. RTI remains selectable with `x_logical_ctrl(3)`, but the test and both pulse
   development cases currently set it to `.false.`.
9. The case-local opacity hook corrects only the ideal-EOS contribution to
   `lnfree_e` above `logT = 3.7` and smoothly continues the opacity below the
   table edge at `logT = 2.7001`, joining the ordinary table by `2.8001`.
10. The shared inlist keeps outer-2%-by-mass drag available while `v_flag` is
    active and uses artificial viscosity in both phases. Both phase inlists
    enable the radiation-pressure floor for the momentum boundary and select
    the `dPrad/dm` temperature-gradient equation.

The single-star development case and the binary-product pulse case are open
ended (`max_model_number = -1`) and do not install test-suite helpers. The
test retains `max_model_number = 15000`, the required termination code, and
the stop 100 days after the first pulse relaxation.

Primary files:

- `star/test_suite/ppisn/inlist_ppisn`
- `star/test_suite/ppisn/inlist_hydro_on`
- `star/test_suite/ppisn/inlist_hydro_off`
- `star/test_suite/ppisn/src/run_star_extras.f90`
- `star/dev_cases_test_TDC/ppisn_models/dev_TDC_single_star_ppisn`
- `star/dev_cases_test_TDC/ppisn_models/dev_TDC_ppisn_from_binary`
- `star/dev_cases_test_TDC/ppisn_models/binary_ppisn_progenitor`

## 7. Split/Merge AMR Surface Mass Floor

Status: committed in checkpoint `16433fecd` and installed successfully.

Split/merge AMR now honors the existing `min_surface_cell_dq` control. For
zone 1 its effective mass floor is

```text
dq_min,1 = max(split_merge_amr_dq_min, min_surface_cell_dq),
```

while interior cells continue to use `split_merge_amr_dq_min`. The zone-1
split test uses `dq_min,1`. Before and after ordinary AMR operations, a
dedicated pass conservatively merges zones 1 and 2 until
`dq(1) >= min_surface_cell_dq`. These mandatory merges bypass the
repeated-remesh guard; ordinary merges retain it. The final cell can exceed
the floor by the mass of the last whole neighbor. This prevents metric zoning
from repeatedly splitting the surface toward the global AMR mass floor
without adding a duplicate AMR-specific control.

The current PPISN hydro experiment sets `min_surface_cell_dq = 1d-14`.

Primary files:

- `star/private/adjust_mesh_split_merge.f90`
- `star/defaults/controls.defaults`
- `star/test_suite/ppisn/inlist_hydro_on`

## 8. Nonnegative-Gas-Pressure Momentum Boundary

Status: committed in checkpoint `16433fecd` and enabled for PPISN hydro.

The optional regular control

```text
floor_momentum_outer_BC_at_Prad = .false.
```

changes only the pressure traction used by `use_momentum_outer_BC`. After the
atmosphere returns its pressure and temperature, the effective momentum
boundary pressure is

```text
Prad_bc = a*T_bc^4/3,
P_momentum = max(P_atm, Prad_bc).
```

The normal atmosphere is therefore unchanged when its implied gas pressure

```text
Pgas_atm = P_atm - Prad_bc
```

is nonnegative. If the hydrostatic atmosphere would require negative gas
pressure, the momentum boundary instead uses `Pgas = 0`. The selected branch
retains its automatic-differentiation derivatives. The atmospheric
temperature equation and `Pextra_factor` are unchanged, and other outer
boundary conditions ignore the control.

For the default Eddington relation with `tau = 2/3` and
`Pextra_factor = 1`,

```text
P_atm - Prad_bc = (tau*g/kap)*(1 - L/Ledd).
```

The floor therefore activates in the super-Eddington regime where the static
atmosphere has no nonnegative-gas-pressure solution. It prevents the missing
pressure from becoming an artificial outward surface-cell force divided by a
small `dm(1)`. It is a boundary closure, not relativistic radiation
hydrodynamics, so ejecta removal and the speed-of-light retry remain separate
issues.

Primary files:

- `star/private/hydro_eqns.f90`
- `star/defaults/controls.defaults`
- `star_data/private/star_controls.inc`
- `star/private/ctrls_io.f90`
- `star/test_suite/ppisn/inlist_hydro_on`

## 9. Conservative `adjust_mass` Velocity Remap

Status: implemented after checkpoint `16433fecd`; compiled and exercised in a
PPISN `u_flag` run through model 4979.

For hydrodynamic mass loss, `do_adjust_mass` now conservatively remaps radial
momentum with MESA's existing quad-precision mass-overlap integrator. The
`u_flag` path uses the ordinary cell masses,

```text
u_new,j = sum_i(overlap_ij*u_old,i)/dm_new,j.
```

The staggered `v_flag` path instead uses its dual control masses:

```text
dmv_1 = dm_1/2,
dmv_k = (dm_(k-1) + dm_k)/2,  k > 1,
v_new,k = sum_i(overlap_ik*v_old,i)/dmv_new,k.
```

The same overlaps transport the old specific kinetic energy. For `u_flag`,
their difference from the kinetic energy of the momentum-projected velocity is

```text
Delta_K_j = dm_new,j*(<u^2/2>_j - <u>_j^2/2) >= 0.
```

For `v_flag`, the analogous deficit is first evaluated on each dual cell and
then assigned to the two adjacent ordinary half-cells. This reproduces the
per-cell kinetic-energy convention in `cell_specific_KE`.

The raw total-energy `eps_mdot` calculation already sees this kinetic-energy
deficit. It is therefore subtracted before thermal leakage and restored once
as a local non-leaking source `Delta_K_j/(dm_new,j*dt)`. The combined source
retains the existing `eps_mdot_factor` scaling. Both `u_flag` and `v_flag`
mass loss are supported. Accretion remains unchanged because its incoming
momentum requires a separately defined boundary state.

Primary files:

- `star/private/adjust_mass.f90`
- `star/private/evolve.f90`
- `star/private/eps_mdot.f90`

Detailed derivation and required run tests are in
`notes/EbF/improve_split_merge_amr_mesh/ppisn_include_kinetic_energy_adjust_mass.md`.

## 10. Shared Minimum `Delta lnR` Mesh Floor

Status: implemented and compiled; not yet tested in an evolutionary run.

Split/merge AMR now honors the existing standard-mesh controls
`mesh_min_dlnR` and `merge_if_dlnR_too_small`. No new controls or restart
state are added.

Before splitting a cell at

```text
r_mid = (r_outer + r_inner)/2,
```

both proposed children must satisfy the standard MESA safety margin

```text
min[ln(r_outer/r_mid), ln(r_mid/r_inner)] >= 2*mesh_min_dlnR.
```

This restriction applies to ordinary metric-selected splits and emergency
`dq` splits. When `merge_if_dlnR_too_small` is true, an initial mandatory
pass merges cells with

```text
ln(r_outer/r_inner) < mesh_min_dlnR.
```

These numerical-validity merges bypass the metric, shock, and repeated-remesh
guards. A cell adjacent to a true center skips the logarithmic floor because
`R_center = 0` has no finite logarithm.

The PPISN hydro inlist enables the repair with the existing default floor:

```text
mesh_min_dlnR = 1d-9
merge_if_dlnR_too_small = .true.
```

This prevents the combined `logR`/`logTau` metric from retaining
photospheric cells whose radii differ by only a few floating-point increments.

Primary files:

- `star/private/adjust_mesh_split_merge.f90`
- `star/defaults/controls.defaults`
- `star/test_suite/ppisn/inlist_hydro_on`

The run diagnosis and profile measurements are in
`notes/EbF/improve_split_merge_amr_mesh/ppisn_audit.md`.

The first attempted continuation did not test this implementation. The PPISN
executable was built at 2026-08-14 19:38, while
`adjust_mesh_split_merge.f90` was modified at 22:47, and the executable does
not contain the `dlnR floor merge` diagnostic string. The continued run still
used the old split/merge AMR code. A fresh build followed by a restart from
model 3200 is required before evaluating this change.

The correct checkout was installed on 2026-08-15 and the PPISN executable was
rebuilt against it. Both `build/star/lib/libstar.a` and
`star/test_suite/ppisn/build/bin/star` now contain the `dlnR floor merge`
path. The remaining validation is the restart from model 3200.

## 11. Conservative Pressure Reconstruction for `u_flag` Splits

Status: implemented behind the PPISN opt-in. The MLT turbulent-pressure
extension has not been compiled or run.

The default-off development control

```text
split_merge_amr_reconstruct_pressure_for_u_flag = .false.
```

enables a bounded pressure correction after the existing conservative split
and final potential-energy radius correction. PPISN enables it for testing.
The target preserves the pre-split limited supported total-pressure gradient:

```text
Delta_P_target = P_inner - P_outer
   = -0.5*grad(P_eos + P_MLT,turb)*Delta_r_parent.
```

One internal-energy transfer `delta` adjusts the outer child while deriving
the inner child from exact conservation:

```text
e_outer = e_outer,0 + delta,
e_inner = (E_internal,total - dm_outer*e_outer)/dm_inner.
```

A bracketed EOS solve chooses `delta`. Both energies are restricted to the
original three-cell energy range, and accepted child pressures and
temperatures must remain inside their original stencil extrema. Invalid EOS
states, an unbracketed target, or failed convergence retain the existing
conservative split. The correction preserves the existing mass, composition,
momentum, and internal-plus-kinetic-energy invariants.

The pressure mismatch includes MLT turbulent pressure evaluated from the
remapped face `mlt_vc` and geometry-derived face density, while the conservative
energy transfer changes only `P_eos`. TDC eddy viscosity remains a momentum
source and is not included as pressure. Split/merge AMR does not currently
support RSP2. Detailed derivation and required tests are in
`notes/EbF/improve_split_merge_amr_mesh/fix_amr_splitting_1.md`.

Primary files:

- `star/private/adjust_mesh_split_merge.f90`
- `star/defaults/controls_dev.defaults`
- `star_data/private/star_controls_dev.inc`
- `star/private/ctrls_io.f90`
- `star/test_suite/ppisn/inlist_ppisn`

## 12. Riemann-Shock Chemical Mixing Barrier

Status: implemented in the current worktree; not compiled or run.

The default-off development controls

```text
Riemann_shock_D_mix_reduction_on = 0d0,
Riemann_shock_D_mix_reduction_full_on = -1d0
```

smoothly reduces only the final symmetric chemical diffusion coefficient at a
compressive `u_flag` Riemann shock. At face `k`, with inner state `L = k` and
outer state `R = k-1`, the detector is

```text
C_k = max(0, (u_L-u_R)/cs_face),
J_k = max(0, P_face/min(P_L,P_R) - 1),
S_k = min(C_k, J_k).
```

`P_L` and `P_R` use the same gravity correction as the HLLC momentum solve,
and `P_face` is recomputed with its scalar contact-pressure formula. This
avoids relying on stale face state after a load or remesh. For onset `S_on`
and positive full-on strength `S_full`, the reduction is the bounded quintic
smootherstep

```text
f_k = 1,                                      S_k <= S_on,
x_k = (S_k-S_on)/(S_full-S_on),               S_on < S_k < S_full,
R_k = x_k^3*(10 - 15*x_k + 6*x_k^2),
f_k = 1 - R_k,
f_k = 0,                                      S_k >= S_full.
```

When enabled, control validation requires `0 <= S_on < S_full`. The exact
onset prevents ordinary faces below both thresholds from receiving even a
small reduction.

The final `D_mix`, `D_mix_non_rotation`, optional `D_mix_rotation`, and `cdc`
are multiplied by `f_k` after all mixing contributions are assembled. The
reduction occurs after deriving `conv_vel` and does not change `mlt_D`,
`mlt_vc`, `mixing_type`, `gradT`, `L_conv`, TDC state, turbulent pressure, or
eddy viscosity. This prevents symmetric abundance diffusion through a strong
shock without inserting a one-face hole in the TDC stress or heat flux.

The PPISN hydro-on inlist enables

```text
Riemann_shock_D_mix_reduction_on = 0.05d0,
Riemann_shock_D_mix_reduction_full_on = 0.1d0
```

and disables the two broad velocity-based full-convection shutdowns by setting
their thresholds to `1d99`. No optional convection or `L_conv` shutdown is
provided because it would defeat the separation this control establishes.

Four standard profile columns expose the calculation:

1. `Riemann_shock_compression`
2. `Riemann_shock_pressure_jump`
3. `Riemann_shock_strength`
4. `Riemann_shock_D_mix_factor`

The PPISN profile list enables all four.

Primary files:

- `star/private/hydro_riemann.f90`
- `star/private/mix_info.f90`
- `star/private/star_profile_def.f90`
- `star/private/profile_getval.f90`
- `star/defaults/controls_dev.defaults`
- `star/defaults/profile_columns.list`
- `star_data/private/star_controls_dev.inc`
- `star/private/ctrls_io.f90`
- `star/test_suite/ppisn/inlist_hydro_on`
- `star/test_suite/ppisn/profile_columns.list`

The derivation, physical scope, and diagnostic matrix are in
`notes/EbF/improve_split_merge_amr_mesh/tdc_shock_convection_limiter_ideas.md`.

## 13. Pole-Free Convective Radiative-Luminosity Split

Status: active current-worktree behavior.

`star/private/hydro_temperature.f90:do1_alt_dlnT_dm_eqn` now retains the
MLT/TDC luminosity decomposition for either sign of `gradr`. The face
radiative gradient is assembled as

```text
gradr = f_gradr*P*kappa*L/(16*pi*c*G*m*P_rad),
L0 = 16*pi*c*G*m*P_rad/(f_gradr*P*kappa),
L_rad = L0*gradT,
L_conv = L - L_rad.
```

For nonzero `gradr`, this is algebraically identical to
`L_rad = L*gradT/gradr`. Evaluating `L0` directly removes the removable
zero-over-zero singularity when `L` and `gradr` pass through zero. It uses the
same interpolated or reconstructed face `T`, `P`, and opacity and the same
rotation and user `gradr_factor` corrections as the MLT/TDC solve, retaining
their automatic-differentiation derivatives. Radiative faces continue to use
`L_rad = L`, and no control or stored state is added. The `L_conv` diagnostic
uses the same pole-free coefficient, so reported radiative and convective
luminosities remain consistent with the temperature equation at a zero
crossing.

The luminosity split tests `mlt_mixing_type`, not the final `mixing_type`.
RTI can replace the latter with `rayleigh_taylor_mixing` when its chemical
diffusion coefficient dominates, even though the MLT/TDC state and `L_conv`
remain active. Using the MLT label keeps the temperature equation consistent
with the closure that supplied `gradT`; RTI internal-energy diffusion remains
a separate flux-divergence source in the energy equation.

TDC can retain a finite convection speed and inward enthalpy flux after the
total luminosity and `gradr` become negative. Ordinary MLT can also remain
convective when a still more negative Ledoux reference gradient satisfies
`gradr > gradL`. The old branch instead set `L_rad = L` for all negative
`gradr`, contrary to the existing `L_conv` calculation.

In the PPISN model-4704 failure, face 530 had `gradr = -179.70`,
`gradT = 0.31158`, and `L = -3.9363d4 Lsun`. TDC gives
`L_rad = +68.25 Lsun` and `L_conv = -3.9432d4 Lsun`, whereas the old
temperature equation used the full negative luminosity as radiative and drove
a one-zone thermal collapse. The opacity used in `L0` is the original MLT/TDC
face opacity. The separate `min_kap_for_dPrad_dm_eqn` floor remains applied to
the `dPrad/dm` transport coefficient exactly as before. That floor and the
radiative flux limiter are not involved in this failure: the face opacity is
about `36 cm^2/g`, and the hypothetical flux-limiter factor is one to machine
precision. The old commented suggestion to multiply the residual by
`gradr_factor` was removed: the effective factor is already included in `L0`,
and applying it again would not preserve `L_rad = L*gradT/gradr`.

The detailed run chronology and equations are in
`notes/EbF/improve_split_merge_amr_mesh/tdc_shock_convection_limiter_ideas.md`.

## 14. Lagged-Midpoint TDC Eddy Viscosity for `u_flag`

Status: active current-worktree behavior.

The direct `u_flag` TDC source now uses one start-of-step cell radius,

```math
r^\star_k = r_{\mathrm{mid,start},k}, \qquad
z_k = \frac{u_k}{r^\star_k},
```

in the eddy-viscous strain, time-centered viscous heating, and momentum-force
denominator. With

```math
\Chi_k = C_k(z_{k-1}-z_k), \qquad
F_k = \frac{4\pi(\Chi_k-\Chi_{k+1})}{r^\star_k},
```

the discrete momentum work telescopes to

```math
\sum_k u_k F_k = -4\pi\sum_k C_k(z_{k-1}-z_k)^2 \le 0
```

apart from boundary terms. `rmid_start` is fixed through the nonlinear solve
and is regenerated for split and merged cells from the final remeshed geometry.
The velocity derivatives remain implicit and block tridiagonal. No control or
new stored state is added, and `v_flag` is unchanged. The implementation is in
`star/private/tdc_hydro.f90`; `p_d_u_div_rmid_start` already provides the
matching profile diagnostic.

## 15. PPISN Rayleigh-Taylor Mixing

Status: active current-worktree behavior.

The PPISN test can enable MESA's Duffell RTI model only while cell-centered
Riemann hydro is active. The local test-case switch is

```text
x_logical_ctrl(3) = .true.  ! use_RTI_during_hydro
```

Setting it false retains TDC and its independently configured eddy viscosity
without adding RTI. `run_star_extras` owns the state transition because
`inlist_hydro_on` is applied through `star_read_controls` and cannot change
the solver variable set:

```text
hydro off:  RTI_flag = false
hydro on:   set u_flag, then set RTI_flag from the PPISN control
new pulse:  enabled RTI initializes alpha_RTI to zero
restart:    enabled RTI is a no-op when already active, preserving alpha_RTI
cut/relax:  disable RTI before disabling u_flag
```

The general `relax_to_star_cut` path now restores an incoming `RTI_flag` after
rebuilding the retained model. Its reconstructed model receives a zeroed
`alpha_RTI` field rather than silently leaving RTI disabled.

The optional behind-shock RTI ramp is now skipped when its configured mass
width is zero. For a positive width, its multiplier is bounded below by zero,
avoiding both the previous division by zero at the default setting and a
negative diffusion coefficient ahead of the detected shock.

Split/merge AMR now remaps `alpha_RTI` with the child masses. For parent value
`a` and child difference `delta_a`,

```math
a_R = a + \frac{dm_L}{dm}\,\delta_a, \qquad
a_L = a - \frac{dm_R}{dm}\,\delta_a,
```

so

```math
dm_R a_R + dm_L a_L = dm\,a.
```

The slope is reduced when necessary to keep both children nonnegative and
inside the original three-cell `alpha_RTI` range. There is no fixed upper
clip; in particular, `RTI_max_alpha` remains a source shutoff scale rather
than a remap bound. The existing merge was already mass weighted.

`inlist_hydro_on` records the MESA-default RTI coefficients and the outward
supersonic source gate. It disables the default absolute-mass RTI boost, which
is tied to the inner 4--5 `Msun` and is not appropriate as an implicit PPISN
calibration. PPISN profiles now include `alpha_RTI`, its source and diffusion
coefficients, RTI momentum and energy terms, the face-velocity kick, RTI
chemical diffusion, and the local Riemann shock diagnostics. A separate
four-panel PGSTAR window plots the RTI growth and decay sources, state and
diffusion coefficients, face-velocity kick, and momentum and energy terms
against zone number. History also records `shock_q`, shock radius, and shock
Mach number.

The one-global-shock cutoff and signed `v/cs` source gate are retained. They
match the original outward CCSN use case but may be ambiguous with detached
PPISN shells. Their alternatives and the required TDC/energy comparisons are
documented in `notes/EbF/improve_split_merge_amr_mesh/ppisn_rti.md`.

Primary files:

- `star/private/remove_shells.f90`
- `star/private/adjust_mesh_split_merge.f90`
- `star/private/mix_info.f90`
- `star/test_suite/ppisn/src/run_star_extras.f90`
- `star/test_suite/ppisn/inlist_hydro_on`
- `star/test_suite/ppisn/profile_columns.list`
- `star/test_suite/ppisn/history_columns.list`
- `docs/source/changelog.rst`

## 16. Compatible Eddy-Viscous Work in the Total-Energy Equation

Status: active branch behavior.

The internal-energy and conservative total-energy forms require different
eddy-viscous source terms. If the energy residual does not contain radial
kinetic energy, as in classic RSP and the time-centered
`P d(1/rho)` formulation used by RSP2 pulsation models, the local source is

```math
S_{q,k}=E_{q,k}.
```

When the `dedt` residual explicitly contains `dKE/dt`, the mechanical work of
the viscous acceleration must also appear. For cell-centered `u_flag`,

```math
S_{q,k}=E_{q,k}+\bar u_k U_{q,k},
\qquad
\bar u_k=\frac{u_k^{n+1}+u_k^n}{2}.
```

For staggered `v_flag`, MESA's cell kinetic energy assigns half of the cell to
each bounding velocity face, giving

```math
S_{q,k}=E_{q,k}
+\frac{1}{2}\left(\bar v_k U_{q,k}
+\bar v_{k+1}U_{q,k+1}\right).
```

The arithmetic velocity midpoint is required because
`(v_new^2-v_start^2)/2 = v_bar*(v_new-v_start)`. The `v_flag` source uses the
same half-cell and mass-correction factors as `cell_specific_KE`. At the
center, the unstored center-face contribution is zero. TDC supports both
velocity formulations; RSP2 receives the corresponding `v_flag` term. No new
control is required because the correction follows the selected energy
equation.

The `v_flag` source value includes both bounding accelerations. Its inner-face
AD expression is shifted with the same block-tridiagonal projection already
used by the face-averaged TDC `Eq`; the unavailable second-neighbor derivative
is not added to the matrix bandwidth.

Implementation:

- `star/private/hydro_energy.f90` determines whether the active residual
  contains `dKE/dt` and adds the matching mechanical work only in that case.
- `star/private/tdc_hydro.f90` remains the authoritative source of the TDC
  accelerations used by momentum and energy.
- `star/private/hydro_rsp2.f90` remains the authoritative source of the RSP2
  face acceleration.

This closes the local accounting gap in which a cell accelerated by eddy
viscosity could pay for its kinetic-energy increase by cooling its internal
energy. The time-centered PdV path is unchanged because adding mechanical work
there would count radial kinetic energy twice.

## 17. Compatible RTI Momentum-Diffusion Energy

Status: active branch behavior.

The PPISN model-4790 collapse is distinct from the TDC eddy-viscous failure.
The failing zone is RTI-mixed but MLT-radiative, with `mlt_vc = 0`. At model
4760, its RTI acceleration requires about `2.60d11 erg/g/s` of kinetic-energy
transfer while RTI internal-energy diffusion supplies only `3.10d10 erg/g/s`.
The measured internal-energy loss is `2.12d11 erg/g/s`.

For cell-centered Riemann hydrodynamics,

```math
F_k=\sigma_k(u_{k-1}-u_k)-\sigma_{k+1}(u_k-u_{k+1})
```

is now shared by the momentum and energy implementations. The total-energy
source includes its exact midpoint work and the symmetric interface
dissipation,

```math
S_k=\frac{\bar u_kF_k+Q_k}{dm_k},
\qquad
\sum_k \left(\bar u_kF_k+Q_k\right)=0.
```

The work term prevents an accelerated cell from paying for RTI-smoothed
kinetic energy with internal energy. The dissipation term converts the lost
resolved kinetic energy to local heat without changing the global energy.
The existing RTI internal-energy diffusion remains separate and unchanged.

Implementation:

- `star/private/hydro_riemann.f90` evaluates the shared AD force and
  interface dissipation.
- `star/private/hydro_energy.f90` adds the compatible source only for the
  cell-centered RTI momentum-diffusion path.

The detailed budget and discrete equations are in `ppisn_rti.md`.

## 18. Positive Optical Depth During Surface Removal and AMR

Status: active current-worktree behavior.

Surface removal previously replaced the surface optical-depth factor with

```math
f_{\tau,\mathrm{new}}
 = f_{\tau,\mathrm{old}}\frac{\tau_{\mathrm{eff}}}{\tau_1}
 = \frac{\tau_{\mathrm{eff}}}{\tau_{\mathrm{base}}},
```

whenever the momentum outer boundary was active.  The local estimate

```math
\tau_{\mathrm{eff}} = \kappa\left[
   \frac{P_{\mathrm{eos}}}{g}
   - f_{P,\mathrm{extra}}\frac{L/M}{6\pi cG}
\right]
```

assumes that the retained cell can support a hydrostatic atmosphere.  It can
be nonpositive in super-Eddington ejecta.  Surface removal now installs the
estimate only when it is finite and positive; otherwise it retains the
incoming `tau_factor`.

The metric-zoning `logtau` component also requires every optical-depth
boundary and `tau_center` to be finite and positive.  If that invariant is
not satisfied, only the `logtau` weight is set to zero for the current remesh
call.  The independently normalized `logR` component remains active with the
same `split_merge_amr_nz_baseline`.  This supplies a deterministic fallback
for old photos or external model states without introducing an arbitrary
optical-depth floor.

The PPISN failure began at model 3921 immediately after a direct surface cut.
The cut changed the radius from about `1.99d4` to `3.53d3 Rsun` and left
`tau(1) = -9.61d-27` while `tau(2) = 1.77d-4`.  The raw logarithm in metric
zoning then made the normalization invalid, after which AMR repeatedly added
up to its configured 1000-zone iteration limit.  The model grew from 1632
zones at model 3921 to 64449 zones at model 4073 even though the saved metric
classified approximately 64200 cells as undersized.

Primary files:

- `star/private/remove_shells.f90`
- `star/private/adjust_mesh_split_merge.f90`

## 19. Photospheric Luminosity and Effective Temperature

Status: active current-worktree behavior.

When `tau_factor < 1`, `get_phot_info` locates the photosphere at
`tau_phot = tau_base` inside the model.  For the cell spanning that optical
depth, it interpolates

```math
f = \frac{\tau_{\mathrm{phot}}-\tau_k}{\Delta\tau_k},
\qquad
L_{\mathrm{phot}} = L_k + f(L_{k+1}-L_k).
```

The interpolated luminosity was then discarded and replaced by
`max(1d0, L(1))`.  This made `photosphere_L` and the blackbody effective
temperature follow the surface luminosity even when the photosphere was far
inside a dynamic model.  A negative surface luminosity therefore produced a
one-erg-per-second photospheric luminosity and a nearly zero reported `Teff`,
even when the luminosity at `tau = tau_base` was positive.

The positive floor now applies to the interpolated luminosity itself:

```math
L_{\mathrm{phot}} \leftarrow \max(1\ \mathrm{erg\ s^{-1}},
L_{\mathrm{phot}}),
\qquad
T_{\mathrm{eff}} = \left(
\frac{L_{\mathrm{phot}}}{4\pi R_{\mathrm{phot}}^2\sigma}
\right)^{1/4}.
```

The surface fallback remains unchanged when `tau_factor >= 1` or no interior
photosphere is found.  The standard `log_L` history column also remains the
surface luminosity `L(1)`.  This is a reporting correction and does not alter
the structure equations, atmosphere boundary, or timestep.

Primary file:

- `star/private/star_utils.f90`

## 20. Reverted or Deferred Work

Status: not active branch behavior.

1. The global piecewise-monotonic cubic remap of `mlt_vc` was committed and
   then reverted. Surviving split/merge faces retain their values; only a new
   split face receives the existing local linear initialization.
2. Raw `logT` metric zoning was removed after it repeatedly refined a persistent
   PPISN shock and collapsed the CFL timestep.
3. Capped `logT` and `logRho` feature metrics remain deferred.
4. Piecewise-linear Riemann face reconstruction remains a design note. The
   bounded pressure-child reconstruction is implemented as described above.
5. The PPISN retained-surface settling, age restoration, conservative
   composition remap, transactional cleanup, and physical outflow-boundary
   findings remain deferred. Only cut indexing, energy diagnostics, and the
   computational extended-layer prune are implemented in the test case.
6. A TDC heat-flux shock limiter remains deferred. The implemented
   composition-only barrier deliberately leaves `L_conv`, `mlt_vc`, and eddy
   viscosity unchanged. Any future heat-flux response requires evidence of
   upstream leakage and a formulation that preserves the monotone local TDC
   solve.
7. The optional Cheng-Shu general-EOS HLLC wave-speed bound was removed. It
   was inactive in the accepted PPISN state and did not prevent the cooling
   failure; the measured defect was missing RTI momentum-diffusion work.

## 21. Verification State

The committed branch previously completed a full MESA installation. For the
changes after checkpoint `16433fecd`:

1. `git diff --check` passes.
2. `fortitude check star/private/tdc_hydro.f90` passes.
3. `sphinx-lint docs/source/changelog.rst` passes.
4. `fortitude check star/test_suite/ppisn/src/run_star_extras.f90` passes.
5. `fortitude check star/private/adjust_mesh_split_merge.f90` passes.
6. A full `./install` completed successfully after the cyclic AMR surface
   mass-floor change.
7. `fortitude check star/private/hydro_eqns.f90 star/private/ctrls_io.f90`
   passes after adding the momentum-boundary pressure floor.
8. A full `./install` completed successfully after adding
   `floor_momentum_outer_BC_at_Prad`.
9. A PPISN smooth-removal run reached model 4979 without the original immediate
   thermal collapse. It removed `1.795 Msun`, then stopped under the previous
   sound-speed gate because the surface Mach number fell to `0.974` even though
   `u/vesc = 6.49`. The sound-speed gate has now been removed, and the PPISN
   inlist instead requires `r > 1d4 Rsun`, `u > 4*vesc`, and positive specific
   total energy. Smooth removal was subsequently deleted; the maintained cases
   use direct removal only.
10. A subsequent direct-removal run shows repeated small cuts rather than one
    clean truncation. From models 4288 through 4363, `star_mass` falls from
    about `51.631 Msun` to `50.505 Msun` while `star_mdot` remains zero and
    `log_R` repeatedly returns to about 4. The cut indexing is consistent; the
    behavior follows from using `R_limit = 1d4 Rsun` as the per-cell stopping
    radius while the newly exposed surface still moves outward at thousands of
    km/s.
11. `git diff --check` passes for the conservative `adjust_mass` remap.
12. A full `./install` completed successfully with the conservative
    `adjust_mass` remap for both velocity formulations; no model run was made.
13. Standalone overlap checks for both `u_flag` and `v_flag` preserve momentum
    to machine precision, leave uniform velocity unchanged, keep remapped
    velocities bounded, and reproduce the per-cell projection-energy deficit.
14. `git diff --check` and
    `fortitude check star/private/adjust_mesh_split_merge.f90` pass after
    adding the shared minimum-`Delta lnR` behavior.
15. A profile-4600 reconstruction check activates the thermodynamic split
    limiter for 437 candidate cells. Momentum and internal-plus-kinetic energy
    close to `2.3d-16` and `4.6d-16`, respectively, and the child internal
    energies remain above the stencil minimum to roundoff.
16. A full `./install` completed successfully with the shared minimum-`Delta
    lnR` behavior and thermodynamic split limiter. No PPISN model was run. The
    PPISN `rn1` script calls `make run`, which will relink its older local
    executable against the rebuilt MESA library before the continuation test.
17. A continuation through model 4867 completed the staged composition
    relaxation but failed every entropy-relaxation attempt with
    `get_T_tau -- L <= 0`, down to `dt < 1d-20 yr`. The first fix was compiled,
    but its temporary `fixed_Tsurf` assignment was overwritten by the following
    control reread.
18. Applying `fixed_Tsurf` after that reread kept it active too early. The next
    continuation failed at composition model 1 with `lambda = 0`, repeated EOS
    density and surface-temperature failures, and a 2328-cell replacement
    instead of the previous 848-cell model.
19. Restricting `fixed_Tsurf` to `do_relax_entropy` restored the successful
    848-cell composition relaxation, but entropy then failed at model 1 with
    alternating temperature limits, EOS failures, and density-domain
    violations down to `dt < 1d-20 yr`. The temporary atmosphere changes were
    reverted; the current code uses the normal atmosphere throughout
    `star_relax_to_star_cut`.
20. A PPISN experiment increased the entropy-relaxation timescale from
    `1d-15 yr` to `1d-9 yr`, with `max_dt_for_relax_entropy = 5d-11 yr`. This
    reduced the instantaneous artificial source by `1d6` while preserving
    `max_dt/timescale = 0.05`. A later implicit-source test used `1d-20 yr`
    and `5d-22 yr`, again preserving the same ratio; these are the current
    PPISN inlist values.
21. The longer entropy timescale completed all 30 requested timescales, but its
    unweighted error plateaued at `0.05687`. The residual is concentrated in
    an over-resolved, low-mass outer layer with sharp retained entropy jumps;
    its approximate mass-weighted error is `1.2d-5`. Ordinary evolution then
    failed at the same layer. More elapsed timescales at fixed relaxation
    strength are not expected to help; staged tightening remains to be
    implemented and tested.
22. Entropy relaxation now supports
    `relax_entropy_use_implicit_source`. It is false by default and enabled in
    the PPISN inlist. The implicit path preserves the existing heating value
    while including its local EOS density and temperature derivatives in the
    structure Jacobian; both energy hooks are restored after relaxation.
23. Entropy relaxation now also supports an opt-in normalized source,

    ```math
    \epsilon_{\rm relax}=\frac{T(S_{\rm target}-S)}{\tau_{\rm relax}}.
    ```

    With fixed composition, `T dS/dt = epsilon_relax` gives

    ```math
    \frac{dS}{dt}=\frac{S_{\rm target}-S}{\tau_{\rm relax}},
    ```

    so the requested timescale is the local entropy e-folding time rather than
    being multiplied by `T*S_target/e`. The default remains the original
    `e*(1-S/S_target)/tau` source. PPISN enables the normalized and implicit
    paths together. During entropy relaxation only, it also sets
    `scale_max_correction = 0.1` and `retry_hold = 5`; both values and both
    energy hooks are restored afterward. The implementation is in
    `star/private/relax.f90`, with controls in the existing `star_job` entropy
    block. A full `./install` and the PPISN executable relink completed
    successfully. No evolutionary model has been run with this implementation.
24. For the pressure-child reconstruction, `git diff --check` and
    `fortitude check star/private/adjust_mesh_split_merge.f90
    star/private/ctrls_io.f90` pass. No MESA compilation or model run has been
    performed for this current-worktree change.
25. For the Riemann-shock chemical mixing barrier, `git diff --check`, targeted
    `fortitude`, and `linters/check_columns.py` pass. The development-control
    defaults/include/I/O checks are complete; `check_defaults.py` reports only
    the pre-existing regular-control `overshoot_f2` mismatch. No MESA
    compilation or model run has been performed for this change.
26. PPISN profiles through model 4548 show the shock diagnostics confined to
    one outward-moving compressive front, generally over one to five faces.
    Pressure, temperature, entropy, and velocity change coherently there. The
    strongest saved `Riemann_shock_strength` is `0.0873`, so the configured
    `full_on = 0.1` never produces complete suppression; the minimum saved
    diagnostic factor is `0.108`. At model 4548 the weakened front has
    `S = 0.0470` and factor one, confirming the exact onset. The current run
    validates localization, but its already-smooth composition does not test
    preservation of a sharp abundance interface.
27. For the pole-free convective luminosity split, `git diff --check`, targeted
    `fortitude`, and `sphinx-lint docs/source/changelog.rst` pass. A full
    `./install` completed successfully with `MESA_DIR` set explicitly to this
    checkout. No PPISN continuation has been performed for this change.
28. For the lagged-midpoint `u_flag` TDC source, `git diff --check`, targeted
    `fortitude`, and `sphinx-lint docs/source/changelog.rst` pass. A full
    `./install` completed successfully with `MESA_DIR` set explicitly to this
    checkout. No PPISN continuation has been performed for this change.
29. For the RTI restoration, behind-shock ramp, AMR remap, and configurable
    PPISN lifecycle, `git diff --check`, targeted `fortitude`,
    `linters/check_columns.py`, and `sphinx-lint docs/source/changelog.rst`
    pass. The toggle audit covers initial startup, hydro-active and hydro-off
    restarts, transition into Riemann hydro, and transition into cut/relax. No
    MESA compilation or PPISN model run has been performed for these changes.
30. For compatible eddy-viscous work in total-energy equations,
    `git diff --check`, `fortitude check star/private/hydro_energy.f90`, and
    `sphinx-lint docs/source/changelog.rst` pass. A full `./install` completed
    successfully with `MESA_DIR` set explicitly to this checkout. No PPISN
    continuation has been performed.
31. The model-4790 continuation identifies RTI momentum diffusion as the
    remaining one-cell cooling source. Compatible midpoint work and local
    dissipation are implemented for `u_flag`. The shared force reproduces the
    saved zone-406 `dudt_RTI` to roundoff, and the model-wide midpoint-work
    and dissipation sums cancel to `2.7d-16` relative error. `git diff
    --check`, targeted `fortitude`, and `sphinx-lint` pass. No MESA
    compilation or continuation has been performed for this correction.
32. The experimental general-EOS Riemann wave-speed control and implementation
    were removed after they did not prevent the cooling failure. The original
    Davis acoustic bounds remain. `git diff --check`, targeted `fortitude`,
    and `sphinx-lint` pass; no installation was performed after the removal.
33. The old development PPISN directory was replaced by isolated-star,
    binary-product pulse, and binary-progenitor cases. The two pulse
    development sources are byte-identical; the test copy differs only in
    test-suite bookkeeping and its 100-day success stop. Targeted `fortitude`
    and `git diff --check` pass. No compilation or model run was performed for
    this synchronization.
34. The retained hydro-on and hydro-off inlists explicitly set
    `restore_mesh_on_retry = .false.`. The earlier five-step mesh-restoration
    experiment was removed; `split_merge_amr_avoid_repeated_remesh` remains a
    separate within-pass AMR safeguard. The three phase inlists are
    byte-identical and `git diff --check` passes.
35. For the surface-removal optical-depth guard and invalid-`logtau` AMR
    fallback, `git diff --check`, targeted `fortitude`, and `sphinx-lint`
    pass. Applying the fallback to model 4073 disables `logtau` as expected
    and classifies 64212 of 64449 cells as undersized under the remaining
    normalized `logR` metric. A full `./install` completed successfully with
    `MESA_DIR` set explicitly to this checkout. No PPISN continuation was run.
36. For the photospheric-luminosity reporting correction, `git diff --check`,
    `fortitude check star/private/star_utils.f90`, and
    `sphinx-lint docs/source/changelog.rst` pass. No compilation or model run
    was performed while the current PPISN continuation was active.
37. The default-off `constant_L` control is restored for idealized
    hydrodynamic tests. It replaces the temperature-gradient residual by

    ```text
    R_L(k) = [L(k) - L(k+1)]/L_scale,
    L(nz+1) = L_center,
    L_scale(k) = max[1, |L_center|, |L_start(k)|, |L_next,start(k)|],
    L_next,start(k) = L_start(k+1) for k < nz and L_center for k = nz,
    ```

    and uses the same equation at the surface instead of a temperature
    boundary condition. The `delta_lgL` timestep limit is skipped because the
    prescribed zero luminosity can otherwise change sign at roundoff level. It
    does not select a momentum boundary condition or alter model-file storage.
    A full `./install` completed successfully with this timestep correction;
    the benchmark has not been rerun. When an explicit momentum boundary is
    selected, `PT_eqns_surf` skips the unused atmospheric pressure-temperature
    closure; this avoids evaluating `log(T_surf + dT0)` in cold gamma-law
    tests. The corresponding Newton-correction weight is

    ```text
    w_L(k) = 1/max[1 erg/s, 0.1*|L_start(1)| + |L(k)|].
    ```

    The floor prevents zero-luminosity states from forming `0*infinity` in
    `sizeB`; the magnitude prevents inward luminosity from cancelling the
    correction scale. The restored Sedov and Noh cases set
    `include_L_in_correction_limits = .false.`, matching their legacy behavior:
    the prescribed zero luminosity remains in the Newton equations and
    residual tests but is excluded from the correction-size test. The cases use
    test-local public hooks for

    ```text
    P = (gamma - 1)*rho*e,  cgrav(k) = 0,  opacity(k) = 0.2,
    ```

    with `gamma = 1.4` and `1.66667`, respectively. Their legacy r15140 model
    controls also exempted gamma-law models from the hydro matrix temperature
    bounds. Their pressure floors correspond to auxiliary temperatures as low
    as `logT = -16.1`, so the current cases set
    `hydro_mtx_min_allowed_logT = -99` locally. The post-hydrodynamic
    convergence check now uses this existing control instead of imposing a
    separate `logT = 1` floor. Its unchanged default preserves the physical
    temperature floor for normal stellar models. A full `./install` completed
    successfully with this correction; the benchmarks have not been rerun.
    Their legacy controls skipped Brunt evaluation whenever
    `gamma_law_hydro > 0`. The
    test-local replacement sets `calculate_Brunt_B = .false.` and
    `calculate_Brunt_N2 = .false.` because these zero-gravity, uniform-pressure
    states otherwise evaluate

    ```text
    B(k) = [ln P(xa_k) - ln P(xa_{k-1})]/[Delta ln P(k)*chiT(k)]
    ```

    with `Delta ln P = 0`; buoyancy is neither defined nor used with
    `MLT_option = 'none'`. Their legacy r15140 model
    files and build trees were removed; each `rn` now generates a fresh
    current-format model under `mods/`. Inlists use current indexed includes,
    `load_model_filename`, `log_center_temp_upper_limit`, and
    `MLT_option = 'none'`. Their history and profile column files are current
    default catalogs with the time, energy, radius, hydrodynamic state, and
    face-state columns needed for diagnostics and PGSTAR enabled.
    Implementation references are `star/private/hydro_temperature.f90`,
    `star/private/hydro_eqns.f90`, and the two benchmark
    `src/run_star_extras.f90` files. Targeted `fortitude`, `check_pgstar.py`,
    `sphinx-lint`, shell syntax checks, and `git diff --check` pass. A full
    `./install` completed successfully with `MESA_DIR` set explicitly to this
    checkout, and both benchmark work directories compile and link. The Noh
    model now loads and reaches its first solve; that solve exposed the
    zero-luminosity correction-weight guard added afterward. A full
    `./install` completed successfully with the guard, and both benchmarks now
    exclude their prescribed zero luminosity from correction-size limits. The
    Noh model has not been rerun with these changes, and the Sedov model has not
    been run. Both profile-panel plots save 27-inch files every 10 model
    timesteps in their local `png` directories. For its first 1000 accepted
    models, the Noh case also writes `noh_energy_balance.txt` with the
    cell-by-cell discrete balance

    ```text
    Delta U_k = dm_k*(e_k - e_k,start),
    Delta K_k = dt*dm_k*dkedt_k,
    Delta W_k = dt*dm_k*dwork_dm_k,
    B_k = Delta U_k + Delta K_k + Delta W_k.
    ```

    Here `dkedt` uses MESA's `u_flag` kinetic-energy quadrature and
    `dwork_dm` is the divergence of the same `A*P_face*u_face` work used by
    the energy residual. The file includes the unscaled `ergs_error` for a
    direct closure check, per-model sums, and only accepted steps, so retries
    cannot duplicate records. The balance covers the structure solve after
    `energy_start` is established; energy changes from a preceding mesh remap
    remain part of phase-one accounting rather than this diagnostic.
    `plotter/plot_energy_balance.py` streams this potentially large file to
    produce integrated closure histories, central wall-heating histories,
    selected radial balance profiles, and a compact per-model CSV summary.
    The history plots normalize the global and cell-local closure errors by

    ```text
    E_step = |Delta U| + |Delta K| + |Delta W|,
    ```

    and compare the postshock density, pressure, and specific energy with the
    analytic Noh values. In the first 1000 accepted models, ending at
    `t = 0.04654 s`, the maximum global relative balance is `2.14d-9` and the
    maximum cell-local balance is `1.16d-5`. The independent balance differs
    from `ergs_error` by no more than `2.40d-16` of the local step energy.
    Thus the recorded wall-heating structure is not caused by a missing term
    in this discrete energy equation. At model 1000 the central pressure is
    `20.75`, close to the analytic `64/3`, while the central density is
    `208.3` and specific energy is `0.1494` instead of `64` and `0.5`.
    Neighboring postshock material reaches `rho = 43.29`. The error is a
    redistribution of density and internal energy at nearly correct pressure;
    the largest shock-local `|Delta U|` follows the analytic shock at
    `r/r_s = 1.034`. The entropy proxy

    ```text
    K = P/rho^gamma
    ```

    isolates the wall-heating error. At model 1000 the low-density shell has
    `K/K_Noh = 1.866`, while the innermost high-density cell has
    `K/K_Noh = 0.136`. The high-entropy shell remains near
    `r = 8.1d-4 cm`: its position changes from `r/r_s = 0.219` at model 250
    to `0.052` at model 1000 as the shock expands. The excess entropy was
    generated during shock formation at the origin and remains behind there;
    it is not continuing to follow the shock front. This is the Noh
    wall-heating error despite the accurately closed total-energy update.
    The smallest direct numerical remedy is a conservative artificial heat
    flux during strong compression,

    ```text
    L_h = -4*pi*r^2*rho*chi_h*de/dr,
    chi_h = C_h*Delta r*a_signal*S_shock.
    ```

    The same face luminosity must enter the two neighboring energy equations
    with opposite signs and full implicit derivatives. It would redistribute
    the excess internal energy generated near the wall without changing total
    energy or momentum. Increasing artificial viscosity is not equivalent and
    can increase wall heating. Finer central zoning can confine the artifact
    to less mass, but does not remove the underlying shock-startup entropy
    error. A compatible momentum-energy discretization is the fundamental but
    substantially larger alternative.
    The current benchmark resolution is not identical to the MESA IV figures.
    Its published Noh calculation used 10,000 cells and a fixed
    `dt = 5d-6 s`. The current case starts with 1000 cells, varies between
    1000 and 2039 cells, and uses `max_timestep = 5d-5 s`; its final shock cell
    has `Delta r = 6.81d-4 cm`, compared with the published run's nominal
    `1d-4 cm` initial spacing. It is therefore about 6--7 times coarser at the
    shock and takes about one tenth as many timesteps. The current Sedov case
    starts with 1000 cells, remeshes toward a 2000-cell target, and finishes
    with 2242 cells. Its final shock-region median `Delta r` is
    `4.23d-4 cm`, with a minimum of `1.18d-4 cm`. This matches the retained
    MESA IV test's 2000-cell spatial target. The active
    `dt_div_min_dr_div_cs_limit = 0.15` now matches the value marked for the
    paper calculation; the faster `0.75` development value remains commented
    in the inlist.
    Of the branch reconstruction changes, only the conservative `u_flag`
    split/merge velocity remap is active in the present Noh and Sedov inlists
    (`star/private/adjust_mesh_split_merge.f90`). A merge now uses

    ```text
    u_new = (dm_1*u_1 + dm_2*u_2)/(dm_1 + dm_2)
    Delta KE = 0.5*dm_1*dm_2*(u_1 - u_2)^2/(dm_1 + dm_2)
    ```

    and adds `Delta KE` to internal energy. A split reconstructs velocity while
    preserving child momentum and subtracts the newly resolved kinetic energy,
    with a bound that prevents a new stencil internal-energy minimum. This
    removes remesh-created momentum and total-energy errors and can reduce
    remesh noise in both tests. It does not change the piecewise-constant HLLC
    face states in `star/private/hydro_riemann.f90`, increase shock resolution,
    or correct the Noh wall-heating entropy error. The bounded pressure-child
    reconstruction and metric zoning controls are both false in these inlists,
    while RTI, TDC, and shock mixing reductions are inactive. Enabling
    pressure-child reconstruction may reduce split-induced
    pressure noise in Sedov, where the pressure gradient is large, but should
    have little effect on Noh's nearly constant post-shock pressure and is not
    a substitute for the published resolution and timestep.
38. The current-format weak-shock case restores the fixed-luminosity
    black-body surface closure used by the MESA IV calculation. The removed
    legacy controls are replaced by a test-local `other_surface_PT` hook in
    `star/dev_cases_test_TDC/dev_TDC_weak_shocks/src/run_star_extras.f90`:

    ```text
    T_surf^4 = L_fixed/(4*pi*sigma*R_surf^2),
    dlnT_surf/dL = 0,
    dlnT_surf/dlnR = -1/2.
    ```

    `L_fixed` is initialized from the loaded surface luminosity, stored in
    the case's restart data, and reset from the newly exposed surface after
    direct removal. The hook supplies the temperature boundary while the
    existing compression or fixed-velocity momentum boundary remains in
    effect. This prevents the ordinary `T_tau` atmosphere from rejecting
    negative trial luminosities during shock injection without changing
    normal MESA surface boundaries. The Zenodo weak-shock inlists used
    `split_merge_amr_MaxLong = 1.1` and
    `split_merge_amr_MaxShort = 1.1`. That pair was retained during the
    current-format conversion, but its product is less than two and permits a
    newly split cell to be immediately eligible for merging. The prepare,
    injection, and finish inlists now retain the tight `MaxLong = 1.1` split
    threshold and use `MaxShort = 2.5` for valid hysteresis.

## 22. Superadiabatic Turnover-Limiter Photo Restart State

Status: fixed in the current worktree.

Current-format photos store the applied superadiabatic reduction factor
`Gamma_factor`. The turnover limiter advances its reciprocal,

```text
eta = 1/Gamma_factor,
eta_new = eta_old + f_turnover*(eta_inst - eta_old),
```

so the saved `Gamma_factor` is the accepted state required by the next step.
A value of one means no reduction.

The PPISN restart from `photos/x500` exposed an initialization-order bug.
The saved model-11500 photo contains

```text
1 <= superad_reduction_factor <= 89.3174249452.
```

The uninterrupted `x050` through `x500` sequence changes the maximum smoothly
from `89.3174123226` to `89.3174249452`. After restarting from `x500`, the new
model-11550 `x550` instead had a maximum of `1.13735977748`. Binary inspection
of the photo record confirms that `x500` contains the correct factor, so the
state was lost after reading rather than omitted during photo output.

The restart order is:

```text
do_load1_star
  -> read_star_photo
  -> finish_load_model
     -> set_vars
extras_startup
  -> star_read_controls('inlist_hydro_on' or 'inlist_hydro_off')
```

At the time of this failure, both hydro inlists enabled
`use_superad_reduction` and `superad_reduction_use_turnover_limit`, but
`run_star_extras` selected them only after `finish_load_model`. The previous
hold condition also required those controls to be active. During
`finish_load_model`, it therefore reset the restored array to one before
`extras_startup` could restore the phase-dependent controls. The current PPISN
comparison cases place these shared controls in `inlist_ppisn`, but the core
restart fix remains necessary for cases that select them dynamically.

The core fix treats a photo as an authoritative snapshot. When
`have_superad_reduction_factor` is true, `doing_finish_load_model` now holds
the saved factor independently of the current controls. Outside load
reconstruction, the existing condition still holds the factor only while the
turnover limiter is enabled and the accepted state has not yet been copied to
the old arrays. Thus:

- current-format photo state survives load reconstruction;
- `.mod` loading is unchanged because it does not set
  `have_superad_reduction_factor`;
- a run with superadiabatic reduction disabled still resets the factor to one
  on its first real evolution step;
- MLT, TDC, hydro, diffusion, and other restart state are unchanged.

Primary files and call sites:

- `star/private/photo_in.f90:134-143` restores the current and old factors;
- `star/private/init.f90:899-909` reads the photo and calls
  `finish_load_model`;
- `star/private/read_model.f90:161-168` rebuilds variables during load;
- `star/job/run_star_support.f90:459-480` shows that `extras_startup` follows
  `do_load1_star`;
- `star/private/turb_support.f90:330-339` implements the corrected hold;
- the historical PPISN `src/run_star_extras.f90:674-689` selected the
  phase-dependent controls after model loading.

`./install` completed successfully on 2026-08-22, including the `star`,
`astero`, and `binary` checker comparisons. The physical validation still
requires a new user-run restart from `x500`. The first displayed model and the
next photo should retain the pre-restart factor profile up to the physical
turnover update, rather than resetting the active region to one.

## 23. PPISN Ideal-EOS Electron Component in Opacity

Status: implemented in the current worktree; static Fortran lint passes, but
the change has not been compiled or run.

The ideal EOS is the final coverage fallback below FreeEOS, OPAL/SCVH, and
HELM. It sets a nominal electron number density `xnefer = 1d-20`, so its
opacity input is

```text
lnfree_I = ln[1d-20/(N_A*rho)].
```

This neutral placeholder entered the logarithmic EOS blend and collapsed the
Compton opacity in the dilute model-11550 layer. The PPISN `other_kap`
callback now replaces only the weighted ideal contribution above
`logT = 3.7`:

```text
Delta_I = ln(Ye) - lnfree_I,
lnfree_kap = lnfree_e + f_ideal*Delta_I.
```

The existing FreeEOS/HELM-to-ideal blend supplies the transition width. No
additional density boundary, opacity floor, or ionization model is introduced.
The matching derivatives are

```text
dlnfree_kap/dlnrho = dlnfree_e/dlnrho
   + (df_ideal/dlnrho)*Delta_I + f_ideal,
dlnfree_kap/dlnT = dlnfree_e/dlnT
   + (df_ideal/dlnT)*Delta_I.
```

The complete EOS derivative vectors were already retained in
`s%d_eos_dlnd(:,k)` and `s%d_eos_dlnT(:,k)`, including `i_frac_ideal`; no new
`star_info` arrays or restart state were needed. The correction changes only
the local inputs passed to `kap_get` and leaves stored EOS thermodynamics and
`s%lnfree_e` unchanged.

Primary files:

- `star/dev_cases_test_TDC/ppisn_models/dev_TDC_ppisn_from_binary/src/run_star_extras.f90`
- `star/dev_cases_test_TDC/ppisn_models/dev_TDC_single_star_ppisn/src/run_star_extras.f90`
- `star/test_suite/ppisn/src/run_star_extras.f90`
- `notes/other_kap_hot_ideal_free_e_fix.md`

## 24. PPISN Low-Temperature Opacity-Table Continuation

Status: implemented in the current worktree; static Fortran lint passes, but
the change has not been compiled or run.

The late PPISN failure at model 17643 repeatedly evaluated zone 163 at the
opacity table's lower temperature boundary, `logT = 2.7`. Retrying at smaller
timesteps returned to the same boundary. The case-local `my_other_kap_get` now
continues the lowest valid opacity through that numerical table edge without
changing the EOS state or core MESA opacity behavior.

Let

```text
x = log10(T),
x0 = 2.7001,
x1 = 2.8001,
u = (x - x0)/(x1 - x0),
w = 3*u^2 - 2*u^3.
```

With `ell = ln(kap)` and `ell0` evaluated at the same density and composition
but at `x0`, the returned opacity is

```text
ell_ext = ell0,                     x <= x0,
        = (1 - w)*ell0 + w*ell(x),  x0 < x < x1,
        = ell(x),                   x >= x1.
```

Thus the opacity is flat with respect to temperature below `x0`, joins the
normal table with a continuous first temperature derivative between `x0` and
`x1`, and is exactly unchanged above `x1`. The boundary opacity is still
re-evaluated at the current density and composition, so this is not a global
constant opacity floor.

In the transition interval, the callback returns

```text
dell_ext/dlnrho = (1 - w)*dell0/dlnrho + w*dell/dlnrho,

dell_ext/dlnT = w*dell/dlnT
   + [6*u*(1 - u)/(ln(10)*(x1 - x0))]*(ell - ell0).
```

Composition partials and `kap_fracs` use the same linear blend weight. Below
`x0`, `dlnkap/dlnT` is zero while the boundary evaluation's density and
composition partials are retained. An error from the boundary `kap_get` call
is propagated rather than hidden.

The existing unbound `k == 1` electron-scattering override retains precedence.
The continuation has no control or saved state, so fresh starts, photo
restarts, retries, and remeshing all use the current trial-cell state.

Primary files:

- `star/dev_cases_test_TDC/ppisn_models/dev_TDC_ppisn_from_binary/src/run_star_extras.f90`
- `star/dev_cases_test_TDC/ppisn_models/dev_TDC_single_star_ppisn/src/run_star_extras.f90`
- `star/test_suite/ppisn/src/run_star_extras.f90`
- `notes/other_kap_hot_ideal_free_e_fix.md`

## 25. Superadiabatic-Reduction Temperature Cutoff

Status: implemented in the current worktree; not compiled or run.

The regular control

```text
superad_reduction_max_logT = 7d0
```

limits superadiabatic reduction to temperatures below

```text
log10(T_start/K) < superad_reduction_max_logT.
```

For `k > 0`, the default therefore disables the reduction at and above
`T_start = 10^7 K`. `turb_support` uses the fixed start-of-step temperature,
so a face cannot cross the cutoff during Newton iterations. Calls with `k=0`
have no stored start-of-step face and retain the previous behavior. Above the
cutoff, the applied factor is one,

```text
Gamma_factor = 1,
gradr_scaled = gradr,
```

and no reduced-gradient MLT or TDC reevaluation is performed. Photo-load
reconstruction continues to preserve the saved factor until phase-dependent
controls are restored; the first real evolution step resets inactive zones to
one.

Both active PPISN development cases and the PPISN test-suite case set the
control explicitly to `7d0` in their `inlist_ppisn` files.

Primary files:

- `star/defaults/controls.defaults`
- `star_data/private/star_controls.inc`
- `star/private/ctrls_io.f90`
- `star/private/turb_support.f90`
- `docs/source/changelog.rst`
- `star/dev_cases_test_TDC/ppisn_models/dev_TDC_ppisn_from_binary/inlist_ppisn`
- `star/dev_cases_test_TDC/ppisn_models/dev_TDC_single_star_ppisn/inlist_ppisn`
- `star/test_suite/ppisn/inlist_ppisn`

## 26. PPISN Test-Suite Synchronization

Status: implemented in the current worktree; static checks pass, but the case
has not been compiled or run.

The single-star development case is now the source of truth for the PPISN
test-suite physics and controls. The two development copies of
`src/run_star_extras.f90` are byte-identical, and all three pulse cases use
byte-identical `inlist_hydro_on` and `inlist_hydro_off` files. The test extras
differs only by retaining `test_suite_extras`, the 100-day success stop, and
the associated test flag. The obsolete `inlist_after_first_pulse` was removed
because the post-pulse `v_flag` timestep cap is applied directly by
`run_star_extras`.

The test-suite case retains only these test-specific controls:

- `max_model_number = 15000` provides a bounded failure mode;
- `x_logical_ctrl(1) = .true.` stops 100 days after the first pulse
  relaxation;
- `required_termination_code_string` verifies that success condition;
- `profile_interval = 200`; the history interval remains one model for test
  diagnostics.

The success condition is

```text
ix_num_relaxations = 1
and star_age - star_age_at_relax > 100 days.
```

The model-building phase, `rn` TestHub workflow, initial `v_flag` activation,
and saved-model handoff remain unchanged.

`git diff --check`, targeted `fortitude`, and `sphinx-lint` pass. The global
defaults consistency script still exits nonzero for the pre-existing
`overshoot_f2` controls-I/O mismatch; it reports no missing mapping for
`superad_reduction_max_logT`.

Primary files:

- `star/test_suite/ppisn/inlist_ppisn`
- `star/test_suite/ppisn/inlist_pulses`
- `star/test_suite/ppisn/inlist_hydro_on`
- `star/test_suite/ppisn/inlist_hydro_off`
- `star/test_suite/ppisn/src/run_star_extras.f90`
- `star/test_suite/ppisn/README.rst`

## 27. PPISN Rates, Output Cadence, and Progenitor Convection

Status: implemented in the current worktree; static input checks complete.

All four PPISN cases now use the deBoer
`r_c12_ag_o16` special-rate table and otherwise retain MESA's default rates:

```text
num_special_rate_factors = 1
filename_of_special_rate(1) =
   c12ag_deboer_sigma_0p0_2000_Tgrid.dat
```

The CF88 triple-alpha block remains beside it as commented configuration. The
shared binary-progenitor inlist applies the same rate to both stars, and
`basic.net` contains `r_c12_ag_o16`. All four cases use
`profile_interval = 200`; both stars in the binary progenitor set that value
in their per-star inlists.

The binary progenitor now also selects the PPISN convection choices that were
explicitly requested:

- `mixing_length_alpha = 2d0`;
- `MLT_option = "TDC"`;
- `include_mlt_corr_to_TDC = .false.`;
- `TDC_include_eturb_in_energy_equation = .true.`;
- `harmonic_dissipation_length_beta = 1d0`.

The PPISN test, three PPISN dev cases, and the 20M test-suite reference all
select

```text
new_net_name = approx21_cr60_plus_co56.net
adv_net = approx21_cr60_plus_co56.net
mass_fraction_limit_for_Skye = 1d-5.
```

Setting both network names keeps direct network changes and automatic
advanced-network selection consistent. All five TDC configurations include
turbulent energy in the energy equation. The binary progenitor intentionally
retains `thermohaline_coeff = 1d0`; the pulse models retain zero.
Both progenitor stars start with `rotation_flag = .false.`, and tidal
synchronization is disabled consistently with `do_tidal_sync = .false.`.

The three PPISN dev run scripts now select their top-level inlists explicitly
with `MESA_INLIST` and no longer source `star/test_suite/test_suite_helpers`.
The two star cases run their model-building and pulse headers directly; their
local `run_model` functions retain restart cleanup and verify that each
requested output model was created. The binary case runs `inlist_binary`
directly, which now includes `inlist_pgbinary`. The dev-only `rn1` scripts and
the binary case's redundant top-level `inlist` were removed. The dev `re`
scripts use the same explicit inlists, and the obsolete `ck` scripts and
`test_suite_extras` hooks were removed from the dev cases. The PPISN test-suite
workflow continues to use `do_one` and `test_suite_extras` for TestHub
bookkeeping and its 100-day success termination.

The README files for the PPISN test, isolated-star dev case, binary-product
pulse dev case, and binary progenitor now describe their current directory
roles, `MESA_INLIST` run and restart paths, model handoffs, and termination
differences. Their update metadata is `2026-08 by EbF`. In particular, the
documentation records the isolated-star mass as `72.5 Msun`, identifies
`MASS2accretor_final.mod` as the binary-derived pulse input, and documents the
binary progenitor directory that previously had no README. The README text
focuses on the physical calculation and model handoffs rather than TestHub
details or output cadence.

Its helium-burning core also uses the PPISN exponential-overshoot parameters,

```text
(f, f0)_burn_He,core = (0.01, 0.005).
```

MESA selects the first matching overshoot rule for a convective boundary, so
the helium-core rule precedes the progenitor's existing catch-all core rule.
The latter retains `(f, f0) = (0.0415, 0.008)` for non-helium core-burning
phases.

The audit found additional physical differences that remain unchanged pending
a deliberate choice:

- the progenitor uses instantaneous superadiabatic reduction with
  `superad_reduction_diff_grads_limit = 1d-3`, while PPISN uses `1d-4` and
  the turnover limiter;
- thermohaline mixing intentionally remains enabled in the binary progenitor
  and disabled in the pulse models;
- the wind prescriptions differ across the binary handoff, and PPISN adds
  hydro-specific opacity, velocity, and ejecta controls that should not be
  copied into the progenitor without a phase-specific reason.

Primary files:

- `star/test_suite/20M_pre_ms_to_core_collapse/inlist_common`
- `star/test_suite/ppisn/inlist_ppisn`
- `star/dev_cases_test_TDC/ppisn_models/dev_TDC_ppisn_from_binary/inlist_ppisn`
- `star/dev_cases_test_TDC/ppisn_models/dev_TDC_single_star_ppisn/inlist_ppisn`
- `star/dev_cases_test_TDC/ppisn_models/binary_ppisn_progenitor/inlist_both`
- `star/dev_cases_test_TDC/ppisn_models/binary_ppisn_progenitor/inlist_binary`
- `star/dev_cases_test_TDC/ppisn_models/binary_ppisn_progenitor/rn`
- `star/dev_cases_test_TDC/ppisn_models/binary_ppisn_progenitor/re`
- `star/dev_cases_test_TDC/ppisn_models/binary_ppisn_progenitor/inlist1`
- `star/dev_cases_test_TDC/ppisn_models/binary_ppisn_progenitor/inlist2`
- `star/dev_cases_test_TDC/ppisn_models/dev_TDC_ppisn_from_binary/rn`
- `star/dev_cases_test_TDC/ppisn_models/dev_TDC_ppisn_from_binary/re`
- `star/dev_cases_test_TDC/ppisn_models/dev_TDC_ppisn_from_binary/src/run_star_extras.f90`
- `star/dev_cases_test_TDC/ppisn_models/dev_TDC_single_star_ppisn/rn`
- `star/dev_cases_test_TDC/ppisn_models/dev_TDC_single_star_ppisn/re`
- `star/dev_cases_test_TDC/ppisn_models/dev_TDC_single_star_ppisn/src/run_star_extras.f90`

`git diff --check` and shell syntax checks for the binary-progenitor run
scripts pass. No MESA model was compiled or run. The imported progenitor
Fortran sources retain 13 pre-existing, autofixable `fortitude` style findings
(one unnamed `end program` and twelve deprecated relational operators); these
do not affect the inlist physics changes above.
