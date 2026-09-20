# RSP2 independent Y implementation

Branch: `EbF/star_lna`. Base: `622075fbf`.
Started 2026-09-18; updated 2026-09-19.

Development update, 2026-09-19: the user has now requested the optional
three-equation implementation. Its active plan and progress record are in
[rsp2_three_equation_implementation.md](rsp2_three_equation_implementation.md).
The first source edits are provisional and incomplete; further source edits
are paused for the full infrastructure/equation audit recorded there. The
one-equation description and earlier validation history below remain the
baseline; earlier statements that the extension is deferred describe that
earlier stage. The new implementation has not been compiled or run.

Status: RSP2 hydro and RSP2 LNA now use the same signed algebraic PII(Y)
relation. The shared cell mixing-length correction was installed; the current
closure cleanup has not been compiled or installed. The user handles all case
relinking and stellar runs. No commit or push has been made.

This note describes the active equations, state integration, and remaining
validation. Superseded experiments and the previous full development record
are preserved under `output/review/rsp2_no_enthalpy_cap_20260919/notes/`.
The proposed three-equation convection model is deferred until the present
RSP2 model has been checked for robust startup, evolution, and restart behavior.
No additional solver slots or evolution equations are introduced here.

## State and scale heights

`w(k)` remains a **cell** variable in RSP2; its specific turbulent energy is
`w(k)**2`. TDC retains its **face** variable `mlt_vc(k)`, with
`w_face = mlt_vc/sqrt(2/3)`. These placements must not be interchanged.
`k=1` is the surface; increasing k goes inward.

RSP2 now solves signed `Y_face(k)` at faces, alongside independent face `L(k)`
and cell `w(k)`. The former Hp slot is `i_Y`, its equation is `i_rsp2_flux`,
and AD slots 19–21 are `i_Y_m1`, `i_Y_00`, and `i_Y_p1`. No extra solver
variable or AD slot was added.

```math
\mathrm{gradT}_k=\mathrm{gradL}_k+Y_{\mathrm{face},k}.
```

`hydro_rsp2:compute_RSP2_gradT` forms gradL from the face grad_ad and selected
Ledoux composition term; the existing dynamical-gradL option is also retained.
`hydro_vars:set_hydro_vars` bypasses the ordinary MLT inversion for RSP2 so it
cannot overwrite the independent Y. It computes RSP2 state before Brunt and
mixing diagnostics. `turb_info` also guards the independent Y against auxiliary
MLT result storage.

Hp and mixing lengths are derived through the actual public functions in
`tdc_hydro.f90`:

- `get_TDC_Hp_face`: reconstructed or ordinary face scale height.
- `get_TDC_mixing_length_face`: the same ordinary/harmonic length used by TDC.
- `get_TDC_mixing_length_cell`: cell EOS pressure and density with the average
  gravity of the bounding faces, followed by the ordinary, alternate-height,
  or harmonic length law. The harmonic radial limit uses the cell's volume
  midpoint. At a full centre the inner gravity is zero; an excised boundary
  uses cgrav(nz)*M_center/R_center**2.

The face helpers are unchanged by the cell-length correction. RSP2 uses the actual
cell mixing length in damping and cell stress, and the actual face length in
fluxes and face stress. `mix_info` uses this same length in D_mix. Hp is cached
for diagnostics, not solved. `RSP2_assume_HSE` was removed from defaults,
control storage, namelist I/O, and all 11 tracked inlists containing it.

The shared helpers retain local AD derivatives. The cell-local length removes
the extra EOS dependency formerly introduced by averaging two face lengths.
Neighbor shifts still use MESA's existing three-zone AD layout; the separate
k+2 dependence of conservative v-grid viscous work remains unresolved.
This is not a claim of an exact full-band Jacobian. User-run derivative and
convergence checks are required.

## Flux and turbulent-energy residuals

`hydro_rsp2:rsp2_flux_residual` supplies the new independent interior row:

```math
R_{\mathrm{flux},k}=
\frac{Lr_k+Lc_k+Lt_k-L_k}
 {\max(|L_{\mathrm{start},k}|,10^{-3}\max_j|L_{\mathrm{start},j}|,1\ {\rm erg\ s^{-1}})}=0,\qquad k>1.
```

This uses the TDC scale from `turb_support:Get_results`, with absolute values
in the global maximum and a 1 erg/s floor for a zero luminosity profile.
At `solver_iter == 0` the same expression uses current L; during solver
iterations it uses L_start and has zero Jacobian derivatives. It is not scaled
by Y or dL/dY. At k=1 this row is `Y_face(1)=0`. The selected
physical surface luminosity/temperature condition remains in its existing row.
Radiative cells still solve the flux row and signed Y; no positivity rule is
applied to Y. The correction norm includes Y with the default unit correction
weight and native column scale `x_scale=max(1,abs(Y_face_start))`, so it measures
`abs(delta_Y)/max(1,abs(Y_face_start))`. B is already column-scaled on entry to
`sizeB`; a further reciprocal-Y weight would incorrectly normalize it twice.
No residual evaluation inverts L to obtain Y.

With the existing RSP2 face interpolation weights, `compute_Lrad_coeff`,
`compute_Lc_terms`, and `compute_Lt` evaluate

```math
Lr_k=\mathrm{Lrad\_coeff}_k(\mathrm{gradL}_k+Y_{\mathrm{face},k}),
\qquad
\mathrm{Lrad\_coeff}_k=
(4\pi r_k^2)\frac{4acT_{\mathrm{face},k}^4}
 {3\kappa_{\mathrm{face},k}\rho_{\mathrm{face},k}Hp_{\mathrm{face},k}},
```

```math
Lc_k=(4\pi r_k^2)\,\overline{T\rho}_k\,w_{\mathrm{face},k}\,PII_k,
\qquad
Lt_k=-\mathrm{RSP2\_alfat}\,(4\pi r_k^2)^2
\overline{\rho^2}_k\Lambda_{\mathrm{face},k}w_{\mathrm{face},k}
\frac{w_{k-1}^2-w_k^2}{\mathrm{dm\_bar}_k}.
```

Here the bars mean the specific existing RSP2 interpolations in those routines:
interpolating T*rho or rho^2 is not replaced by multiplying interpolated
quantities. Face w is interpolated from the two cell w values, not from w^2.
The existing forced-nonturbulent flux and boundary conditions are retained.

`compute_PII_from_Y` evaluates the same signed relation for every Y:

```math
PII_k=\tfrac12\sqrt{2/3}\,
\frac{\Lambda_{\mathrm{face},k}}{Hp_{\mathrm{face},k}}
Cp_{\mathrm{face},k}Y_{\mathrm{face},k}.
```

The same PII_ad enters convective flux and buoyant driving. Its product-rule
variation retains derivatives of Y, Cp, Lambda and Hp and is regular at Y=0:

```math
\delta PII_k=\tfrac12\sqrt{2/3}\left[
\frac{\Lambda_k Cp_k}{Hp_k}\delta Y_k
+Y_k\left(\frac{Cp_k}{Hp_k}\delta\Lambda_k
+\frac{\Lambda_k}{Hp_k}\delta Cp_k
-\frac{\Lambda_k Cp_k}{Hp_k^2}\delta Hp_k\right)\right].
```

Here Cp, Lambda and Hp are the face quantities defined above. Face w is used
in Lc; the cell source uses its own w after averaging face PII/Hp. PII has no
independent time derivative. Making Y a solver variable introduces an
algebraic flux row, not a separate entropy-flux evolution equation.

The cell turbulent-energy residual in `do1_turbulent_energy_eqn` is

```math
R_{w,k}=\frac{w_k^2-w_{k,\mathrm{start}}^2}{\Delta t}
+P_{\mathrm{trb},k}^{\theta}
 \frac{1/\rho_k-1/\rho_{k,\mathrm{start}}}{\Delta t}
+\frac{Lt_k^{\theta}-Lt_{k+1}^{\theta}}{\Delta m_k}
-\mathrm{COUPL}_k-Eq_k=0,
\qquad
\mathrm{COUPL}=\mathrm{Source}-D-Dr.
```

`compute_Source`, `compute_D`, and `compute_Dr` are also used by STAR LNA.
The pressure and luminosity theta weights retain their original controls.
In forced-nonturbulent cells the w row instead enforces w=0. The existing
positive-w startup factorization is retained only when its local factorization
is valid; it is disabled for u_flag with nonzero eddy viscosity because that
heating depends on neighboring cell w. The startup predictor uses the shared
cell mixing length.

## All three temperature rows

`hydro_temperature.f90` is unchanged. Its existing rows consume the new
`gradT_ad` or `Lr_ad`:

1. `do1_gradT_eqn`: `gradT*(lnPeos(k-1)-lnPeos(k)) - (lnT(k-1)-lnT(k)) = 0`.
2. `do1_alt_dlnT_dm_eqn`: radiation-pressure diffusion using RSP2 Lr, including
   the separate optional radiative flux limiter and its existing opacity floor.
3. `do1_dlnT_dm_eqn`: `delm*dlnPdm_qhse*gradT - (T(k-1)-T(k))/Tpoint = 0`.

The optional radiative diffusion factor retains its existing definition.

## Energy forms: keep the four hydro cases separate

The following describes the eddy-viscosity contributions in
`hydro_energy:setup_sources_and_others`. Other physical sources remain in place.
`Eq` always supplies turbulent heating. RSP2 storage is always cell w^2.

| Velocity placement | `use_P_d_1_div_rho_form_of_work` | dedt storage | Work term | Viscous source in energy row |
|---|---|---|---|---|
| face v | true | e + w^2 | local P times volume-rate | Eq |
| face v | false | e + w^2 + kinetic + potential | difference of face P*A*v | Eq + half-face Uq power |
| cell u | true | e + w^2 | local P times volume-rate from Riemann face velocity | Eq |
| cell u | false | e + w^2 + kinetic + potential | difference of Riemann face P*A*u_face | Eq + `(u+u_start)*Uq/2` |

For conservative face work, the v-grid source is exactly the native kinetic
quadrature:

```math
\mathrm{viscous\_work}_k=
\frac{\mathrm{mass\_correction}_k}{2}
\left[\frac{v_k+v_{k,\mathrm{start}}}{2}Uq_k+
\frac{v_{k+1}+v_{k+1,\mathrm{start}}}{2}Uq_{k+1}\right].
```

The mass-correction factor is one when that option is off. The prescribed
inner-boundary acceleration contributes no evolving kinetic energy.
For u_flag, `hydro_rsp2:compute_Uq_dm_cell` returns a **force**, which the Riemann
momentum and energy code divide by dm exactly once. Uq is no longer added to
the Riemann contact velocity.

`eps_grav` is enabled for RSP2 with either velocity grid. Its unscaled energy
residual is

```math
R_{E,k}= -\frac{L_k^{\theta}-L_{k+1}^{\theta}}{\Delta m_k}
+\mathrm{sources}_k+\mathrm{others}_k+\epsilon_{\mathrm{grav},k}
-\frac{w_k^2-w_{k,\mathrm{start}}^2}{\Delta t}
-\mathrm{dwork\_dm}_k=0.
```

For this form `dwork_dm` uses local pressure work **excluding Peos**, because
the EOS compression/internal-energy contribution is already in eps_grav.
There is no explicit kinetic/potential storage or Uq power in this energy row,
even if the PdV control is false. Native eps_grav evaluation, including its
composition, latent-heat, and entropy choices, is reused.

Local work in `eval_simple_PdV_work` uses
`P*([A*v]_k-[A*v]_(k+1))/dm`; the turbulent row uses the density difference.
Their equivalence relies on the converged continuity/volume equation. They
are not identical residuals away from a solution. Likewise, a local PdV row
is not an assertion of exact finite-step total-energy conservation.

## Eq and Uq: placement, time centering, and conservation

The viscosity implementations are separate. Every executable routine in
`tdc_hydro.f90` is unchanged from base `622075fbf`; only the public list exposes
two existing length helpers. RSP2's `compute_Eq_cell`, `compute_Uq_face`, and
`compute_Uq_dm_cell` live in `hydro_rsp2.f90`. They retain the same discrete
power construction with RSP2's own state and forced-zone controls. There is no
shared w adapter and no RSP2 branch inside the TDC viscosity routines.

The distinct placements are:

| Quantity in viscosity | RSP2 | TDC |
|---|---|---|
| cell w for v-grid stress | `wrap_w_00(s,k)` | half the adjacent face `mlt_vc/sqrt(2/3)` values |
| face w for u-grid stress | interpolate cell w with RSP2 face weights | face `mlt_vc/sqrt(2/3)` |
| prescribed inner boundary w for u grid | innermost cell w | existing innermost face `mlt_vc/sqrt(2/3)` |

For v_flag, Chi is cell-centered. With the radius and velocity time choices
implemented in RSP2 `compute_d_v_div_r(...,.true.)`, the discrete identities are

```math
Uq_k=\frac{4\pi(Chi_{k-1}-Chi_k)}
 {r_k^{\mathrm{work}}\Delta m_{\mathrm{face},k}},
\qquad
Eq_k=\frac{4\pi Chi_k}{\Delta m_k}
 \left(\frac{v_k^{\mathrm{work}}}{r_k^{\mathrm{work}}}
       -\frac{v_{k+1}^{\mathrm{work}}}{r_{k+1}^{\mathrm{work}}}\right).
```

The dual face mass contains the existing mass corrections. For u_flag, Chi is
on faces and the cell radius is `rmid_start`, both in strain and force:

```math
\Delta m_k Uq_k=
\frac{4\pi(Chi_k-Chi_{k+1})}{r_{\mathrm{mid,start},k}},
```

```math
Eq_{\mathrm{face},k}=
\frac{4\pi Chi_k}{(\Delta m_{k-1}+\Delta m_k)/2}
\left(\frac{u_{k-1}^{\mathrm{work}}}{r_{\mathrm{mid,start},k-1}}
     -\frac{u_k^{\mathrm{work}}}{r_{\mathrm{mid,start},k}}\right),
\qquad Eq_k=\tfrac12(Eq_{\mathrm{face},k}+Eq_{\mathrm{face},k+1}).
```

The optional inner-boundary face uses half the innermost cell mass and the
prescribed v_center/R_center strain. The outer face stress is zero. Thus the
half-face distribution gives each interior face its full dual mass when
summed over cells, even for unequal dm.

When momentum and heating use the same w, summation by parts gives

```math
\sum_k\Delta m_{\mathrm{velocity},k}
 v_{\mathrm{velocity},k}^{\mathrm{work}}Uq_k
+\sum_k\Delta m_k Eq_k
=-4\pi Chi_{\mathrm{inner}}v_{\mathrm{center}}/R_{\mathrm{center}}.
```

The right side is absent when the boundary stress is zero. This identity uses
v or u and its corresponding mass, not a universal averaging of both grids.
For conservative energy work the velocity is always `(velocity+velocity_start)/2`,
even when optional velocity time centering is off. For local PdV work it is
that average only when velocity time centering is on; otherwise it is the
current velocity. In the latter case the identity is current-velocity power,
not exact finite-step kinetic-energy change. These native choices are retained.

TDC's `TDC_alpha_M_use_explicit_mlt_vc_in_momentum_equation` is retained:
when active it uses old w for momentum and current w for heating. The same-w
cancellation identity does not apply without the corresponding lag difference.
The algebra check tests that difference explicitly rather than hiding it.
RSP2 uses its current cell w in both terms and does not inherit this TDC option.

During review, an overbroad change that removed TDC's forced-zone Uq guards was
found. The combined viscosity implementation was then removed entirely after
the user raised the risk of altering the staggered TDC scheme. TDC retains
all its original routines and guards. RSP2 masks stresses at their proper
locations; the divergence is retained at an active/forced interface. The audit
does not claim an unrestricted TDC conservation identity across every existing
forced-zone choice.

For RSP2 u_flag, `compute_Chi_face` and `compute_Eq_face` accept k=nz+1 for the
prescribed inner boundary; this result uses the last cell's AD slots and must
not be shifted. Interior k+1 results are shifted exactly once by the cell
caller. `compute_d_u_div_r_face` uses rmid_start in both stress and work.
Future consolidation should wait for a separate demonstrated equivalence
of both models' discrete work, boundary, and time choices.

## Remesh, starts, and persistence

- `mesh_adjust:do_Y_face` interpolates signed face Y and sets the surface to zero.
  Existing conservative turbulent-energy remapping operates on dm*w^2.
- `adjust_mesh_split_merge` now supports RSP2. It moves and writes w/Y caches and
  xh slots together; splits/merges conserve dm*w^2. RSP2 composition splitting
  uses a common limiter and mass-weighted child increments, retaining each
  species mass and sum(X). New face Y is a predictor interpolated in mass.
- RSP2 v-grid AMR retains velocity sign and compensates its native half-face
  kinetic-energy change in thermal energy. The u-grid momentum/energy split
  is retained. QHSE pressure reconstruction includes RSP2 turbulent pressure.
  The native potential-energy radius adjustment remains. It preserves AMR's
  midpoint potential-energy definition, which differs from the main face-average
  definition; see the ordinary-remeshing audit below. Exact additional
  mass-correction/rotation energy conservation is not asserted here.
- `tdc_hydro_support:remesh_for_TDC_pulsations`, also used for RSP2 envelopes, conservatively
  remaps w^2 and composition by mass overlap, then intentionally rebuilds the
  thermal state for QHSE including turbulent pressure. This envelope constructor
  is not an energy-conserving thermal remap.
- `Y_face_start`, `w_start`, and `Lt_start` retain the existing accepted/start
  copy path; Hp_face_start was removed. Solver unpacking, retries, and xh copy
  paths now carry signed Y in the former structural slot.
- `.mod` files write/read w and Y_face. Per user instruction there is no legacy
  Hp conversion. Photo format version is 20; the ordinary version check rejects
  old photos. `set_RSP2_flag` preserves an existing u grid and initializes Y
  from gradT-gradL when enabling RSP2.

## STAR LNA

`star_LNA_support` maps the independent Y slot and uses the same flux residual.
Its RSP2 eps_grav inertia includes the EOS entropy coefficients (plus native latent
heat and the PC entropy blend); Peos work is not added again. LNA chooses the
energy form directly from `energy_eqn_option`, including when requested before
the first hydro equation evaluation. `star_LNA_turbulence_closures` reuses the
nonlinear RSP2 Source/D/Dr and uses the same lengths and w placement in its
static stress linearization. At a static zero-velocity background, Eq and
Uq mechanical power are quadratic and have no first-order heating term, while
Uq remains in the linear momentum row. Existing LNA exclusions, including mass
corrections and eps_grav without RSP2, remain.

The manuscript describes the independent-Y formulation, both energy forms,
both velocity grids, and the distinct work choices. Its colored row pairs
retain the model equation and its LNA linearization together. Build/run status
and remesh settings remain in the development notes. The latest 31-page PDF
was rebuilt without LaTeX warnings or unresolved references.
All rendered pages were inspected; the notes and output copies are identical.

## Progress and validation

- [x] Independent Y state, forward flux residual, boundaries and initialization.
- [x] Shared TDC lengths and all three temperature-gradient rows.
- [x] Both energy forms, both velocity grids and both dedt work choices.
- [x] Separate RSP2/TDC viscosity with shared length functions.
- [x] Ordinary remesh, split/merge, envelope construction and restart state.
- [x] Matching RSP2 hydro/LNA PII, Source, D and Dr closures.
- [x] Current manuscript equations and source references revised.
- [x] Non-compiling algebra, source and Fortran lint checks for the current cleanup.
- [x] PDF rebuild, rendered-page overview and detailed review of the revised closure pages.
- [ ] Compile/install this source revision when requested by the user.
- [ ] User-run AD residual checks and stellar startup/evolution/restart tests.

The standalone scripts check algebra and source wiring without executing
Fortran or MESA. Earlier installed revisions passed the package checks;
those results do not validate the current uncompiled changes or establish
stellar runtime stability. Preserve the energy identities below while
checking the unresolved conservative-work Jacobian dependence separately.

User-run validation should cover the following matrix before trusting results:

1. Build the changed `star_data` and `star` libraries in the usual MESA build
   environment, then rebuild the selected RSP2 test work directory.
2. For `dedt`, test both v/u with the PdV control true/false; repeat with velocity
   time centering on/off and eddy viscosity zero/nonzero. Compare energy budgets
   using the correct work identity above, including prescribed boundary power.
3. For `eps_grav`, test v/u and both settings of the PdV control; confirm EOS
   work is not doubled and the row contains no explicit kinetic-energy source.
4. Exercise each temperature row with positive, negative and near-zero Y;
   test damping and Lt on/off.
5. Test zero-w startup, turbulent/radiative transitions, forced outer/inner cells,
   full-center and excised-envelope boundaries, harmonic lengths, and optional
   face reconstruction. Use the existing AD partial-check facilities for Y, w,
   radius, temperature, and density, accounting for the known stencil truncation.
6. Compare uninterrupted evolution with a new-format `.mod` load, photo restart,
   and forced retry. Exercise ordinary remesh, split/merge (including pressure
   reconstruction), and envelope construction; check mass, each species mass,
   dm*w^2, signed Y, positive density/internal energy, and the intended energy budget.
7. Run a TDC control case on each velocity/work combination, including explicit
   momentum w on/off, to confirm the unchanged TDC routines still behave as before in the integrated build.
8. Compare static LNA with the matching nonlinear background for RSP2 v/u and
   both energy forms. Finite-step convergence and mode-growth agreement have
   not been established by source/algebra checks.

## Ordinary remeshing and energy audit, 2026-09-19

The user requested `RSP2_remesh_when_load = .false.` as the repository default;
updated `star/defaults/controls_dev.defaults`. The scheduled shared envelope
remesher remains controlled by `steps_before_remesh_for_TDC_pulsations`.

- [x] Change the legacy immediate-remesh default to false.
- [x] Trace ordinary face-Y interpolation, including sign and boundaries.
- [x] Trace conservative cell w^2 remapping and its reconstruction order.
- [x] Check thermal, kinetic, potential, and turbulent energy accounting.
- [x] Apply ordinary-remesh/integrator fixes and check the scoped diff and numerical invariants.

This work stays in the development notes; the equation manuscript excludes
remeshing settings and build bookkeeping. No MESA build or model run is
authorized by this request.

### Ordinary mesh: corrected cell-energy coordinates

`mesh_adjust:do_mesh_adjust` constructed `xout_old/new` from `dqbar` for the
face-centered velocity remap. `do_u` and `do_etrb` incorrectly passed those
same dual-grid coordinates to their cell-energy overlap calculations while
weighting by physical cell `dq`. Both now pass `old_xq/new_xq` to
`adjust1_u`/`adjust1_etrb`; the unused dual-coordinate arguments were removed.
The overlap roundoff correction now subtracts the excess accumulated mass
instead of adding it. The face-v remap and TDC hydro equations are unchanged.

For RSP2, the physical invariant is

```math
\Delta m_{ji}=\left|\text{new cell }j\cap\text{old cell }i\right|_m,
\qquad
(w_j^{\rm new})^2=
\frac{\sum_i\Delta m_{ji}(w_i^{\rm old})^2}{\Delta m_j^{\rm new}},
\qquad
\sum_j\Delta m_j^{\rm new}(w_j^{\rm new})^2
=\sum_i\Delta m_i^{\rm old}(w_i^{\rm old})^2.
```

Here `w` is a cell value in cm/s and `w^2` is specific turbulent energy in
erg/g, with no factor of one half. The remap integrates piecewise-constant
old cell averages, then takes the nonnegative square root. `u` uses the same
cell intervals for `u^2/2` and retains the donor velocity sign. This ordinary
u remap conserves kinetic energy, not momentum; split/merge AMR instead
preserves cell momentum and transfers resolved/unresolved kinetic energy
to/from thermal energy.

The numerical regression uses old faces `[0, .2, .6, 1]`, new faces
`[0, .1, .2, .6, 1]`, and old `w = [1, 3, 2]`. The old coordinate wiring loses
6.48148% of turbulent energy; the corrected result is `[1, 1, 3, 2]` and
conserves the integral to roundoff. Earlier overlap-algebra checks assumed
physical cell coordinates and therefore did not detect this wiring error.

### Face Y and interpolation order

`mesh_adjust:do_Y_face` already calls `interpolate_vector` with `interp_pm`:
MESA's piecewise-monotonic cubic interpolation, not a mass-conservative remap.
It retains signed values, copies unaffected faces, extends the innermost
stored value to the inner boundary for interpolation, and enforces the
surface residual's `Y_face(1) = 0`. `interp_1d_pm:mk_pmcub` constructs the
limited slopes; `interp_1d_lib:interp_pm` provides the public interface.
This is a shape-preserving cubic, not an unrestricted global spline.

Keep this face interpolation and the conservative w-squared remap for now.
Conservation and higher order are compatible: a future limited linear
reconstruction of cell `w^2` could reduce remap diffusion while preserving
its integral and positivity. Directly applying a point spline to `w` does
not preserve turbulent energy. Raising the reconstruction order should be
motivated by a remesh-resolution study, especially near convective boundaries.

### Energy accounting and its limits

`mesh_adjust:do1_lnT` overlaps EOS thermal energy and, when enabled, corrects
it using `star_utils:cell_specific_KE` and `cell_specific_PE`. With overbars
denoting physical cell-mass overlap averages, its uncapped target is

```math
e_j^{\rm new}=\overline e_j
+\overline{\mathrm{KE}}_j-\mathrm{KE}_j^{\rm new}
+\overline{\mathrm{PE}}_j-\mathrm{PE}_j^{\rm new}.
```

For ordinary Newtonian nonrotating models, `KE = u^2/2` for cell u, and
`KE = (v_k^2 + v_{k+1}^2)/4` for face v, with `v_center` on the inner boundary.
The new w-squared value already equals its old overlap average, so its
contribution cancels separately: it must not be added to EOS thermal energy.
This gives cellwise and integrated conservation of `e + KE + PE + w^2`
when the energy inversion succeeds and its correction is not limited.

Relevant qualifications:

- `mesh_adjust_get_T_from_E` defaults true, but the current
  `star/dev_cases_TDC_Pulsation/dev_TDC_RSP2_Cepheid/inlist_pulses` explicitly
  sets it false. That case uses interpolated lnT; thermal/total energy is
  not constrained by this remap. The user's case setting was not changed.
- `max_rel_delta_IE_for_mesh_total_energy_balance` defaults to 0.05, limiting
  the thermal correction to 5% of the overlap-averaged thermal energy.
  Zero instead requests thermal-energy conservation alone. The existing
  degeneracy guard (`eta_old >= -1d-6`) and failed-EOS-inversion fallback
  also use interpolated lnT. `mesh_adjust_IE_conservation` is not a reliable
  post-EOS total-energy check in those fallback paths.
- The legacy face-v remap preserves its `dqbar*v^2` integral. Its inner
  dual-cell mass differs from the main diagnostic's half-face quadrature.
  `do1_lnT` evaluates the actual new KE, so the uncapped thermal correction
  accounts for this difference; the face-v predictor alone is not proof
  of total-energy conservation.
- Split/merge AMR already conserves RSP2 `dm*w^2`, compensates v-grid KE in
  thermal energy, and transfers u-grid velocity variance consistently.
  New Y faces use a local linear mass-coordinate predictor; retained faces
  are copied. They are not overlap-conserved. Its existing `get_star_PE`,
  `get_cell_energies`, and `get_star_PE_at_fixed_R_center_scale` use
  `-cgrav*m_cell/r_cell`, whereas `star_utils:cell_specific_PE_qp` averages
  `-cgrav*m/r` at the two faces. The AMR radius adjustment therefore preserves
  a different gravitational-energy quadrature. That broader AMR limitation
  remains and prevents an exact total-energy claim against the main diagnostic.
  This audit does not establish relativistic mass-correction or rotational
  energy conservation for AMR.
- `tdc_hydro_support:remesh_for_TDC_pulsations`, the shared envelope remesher,
  already uses physical mass overlap for RSP2 w-squared and monotonic cubic
  interpolation for Y. Its QHSE reconstruction intentionally rebuilds the
  thermal structure and is not a conservative evolution remap.

`star_utils:cell_specific_total_energy`, `eval_total_energy_profile`, and
`eval_deltaM_total_energy_integrals` each include RSP2 `dm*w^2` once, and exclude
the separate TDC turbulent term when RSP2 is active. The latter routine had
a separate missing `sum_dm = sum_dm + dm` after clipping the final cell to
`deltaM`; this is now fixed. Previously a partial-mass integral could keep
adding deeper cells. Existing full-domain callers, including mass-change
accounting, were not affected by that missing stop except for endpoint
roundoff; this is not evidence that their full-star totals were wrong.

Validation: all seven existing groups in `notes/check_rsp2_independent_Y.py`
pass. All four new groups in `notes/check_rsp2_remesh.py` pass: coordinate-bug
regression and mixed split/merge meshes, u/v thermal-plus-turbulent energy
identities, partial-mass energy boundaries, and source wiring. These are
Python numerical/source checks, not execution of the Fortran. The modified
Fortran files parse as Fortran 2008 with includes/preprocessor directives
omitted; `git diff --check` passes. No MESA compilation or model run was done.
The new changes have not been installed; the earlier install predates them.

### Mass integrator: origin and practical scope

The missing accumulation in `star_utils:eval_deltaM_total_energy_integrals`
exists in all of these Git objects: `HEAD` and `origin/EbF/star_lna` at
`622075fbf`, `origin/main` at `fd396fd73`, the saved pre-update branch at
`f6d606939`, and the initial imported source at `13fe30d52`. These are inspected
local Git objects, not a claim about a newly fetched remote tip.
It predates both the local RSP2 edits and the recent branch update.

The fix is warranted: for cell masses `[1, 1, 1]`, unit specific energy, and
requested `deltaM = 1.5`, the former loop integrates 3 instead of 1.5 because
its `sum_dm` remains zero. The restored increment gives 1.5. For the full
mass 3, both give 3. Every current in-repository caller passes `s% mstar`
as the cap; some restrict the cell-index range as well. Thus this is a latent
partial-mass API defect, not an identified energy error in the present pulsation
run. Its current practical impact is small; a genuinely partial-mass caller
could incur a large error. The earlier description tying it to erroneous
mass-change accounting was too broad: current mass-change callers use the full
domain, and the separate `eval_deltaM_total_from_profile` already increments
its mass correctly.

## Y scaling and flux-residual normalization audit, 2026-09-19

The user requested that this audit exclude the negative-w repair issue.
No source, controls, or tolerances were changed, and no build or model run
was performed for this audit.

`solver_support:set_xscale_info` sets the Y column scale to
`max(1d0,abs(Y_face_start(k)))`. The Y correction weight is one and Y is
included in `sizeB`; the formerly proposed extra reciprocal-Y weight is absent.
`star_utils:store_partials` multiplies Jacobian columns by `x_scale`, and
`star_solver:apply_coeff` converts the solved correction back to physical Y.
Thus the measured Y correction is `abs(delta_Y)/max(1,abs(Y_face_start))`.
This changes units for the linear solve and correction criteria; it does not
quantize Y or impose a minimum residual. For absolute Y below one the column
scale is already exactly one.

Separately, `hydro_rsp2:rsp2_flux_residual` uses the interior row

```math
R_k=\frac{L_{r,k}+L_{c,k}+L_{t,k}-L_k}{L_{{\rm scale},k}},\qquad
L_{{\rm scale},k}=\max\left(1,|L_{{\rm start},k}|,
10^{-3}\max_j|L_{{\rm start},j}|\right).
```

The lower bound 1 is in erg/s. The denominator is fixed during Newton
iterations; initialization at `solver_iter == 0` uses current L. The
surface Y row remains `Y_face(1)=0`. The original TDC luminosity scale is in
`turb_support:Get_results`; RSP2 additionally handles signed profiles and
zero luminosity explicitly. No extra flux residual weight is set.

Scaling the equation and its Jacobian by the same fixed nonzero number
leaves the exact Newton step unchanged. It can affect numerical conditioning,
line-search weighting, and acceptance criteria. Multiplying L_scale by ten
immediately changes a printed `1d-7` to `1d-8` for the identical physical
luminosity defect; at fixed gold2 tolerances this relaxes the row's physical
accuracy requirement by ten. Removing the denominator leaves a dimensional
erg/s residual and makes the existing tolerance inappropriate.

The saved `LOGS/profile20.data` has now been overwritten by model 38000 with
350 zones; it is no longer the earlier model-19000 evidence. At the time of
this audit, zones 300--350 have `max(abs(Y)) = 0.389955`, so a nearby step's
native Y scale would be one. The profile does not contain the actual trial
Y_start or L_start. Using its current luminosity as an estimate of L_start,
face 308 has L_scale about 0.33509 Lsun and
`dR/dY ~= 695.35` at fixed other variables, inferred from
`(Lrad/gradT + Lconv/Y)/L_scale` on this negative-Y branch. A residual of
`1d-8` corresponds to a Y correction about `1.44d-11` in that one-variable
estimate. Most sampled inner faces require corrections of order `1d-9`.
Those increments are representable in double precision. Simple flux-sum
roundoff is of order `1d-16` in normalized units; cancellation in gradL+Y
raises the estimate to about `1.2d-13` at face 308, still well below `1d-8`.
These estimates do not bound errors propagated from EOS evaluation, dynamic
gradL, matrix conditioning, or the full coupled update.

The next discriminating diagnostics are the unnormalized luminosity defect,
the fixed L_scale, requested versus realized delta_Y, and the full Jacobian
prediction for the flux row compared with its change after the trial update.
If normalization changes are considered, judge convergence using the original
physical luminosity defect as well as the newly printed residual. No current
stalled-iterate trace was provided for this audit, so the cause of the reported
`1d-7` plateau is not established by the accepted profile.

## Toggleable flux-solver diagnostics, 2026-09-19

- [x] Add default-off `RSP2_report_flux_solver` to controls and namelist I/O.
- [x] Capture the assembled flux Jacobian before factorization and the raw Newton correction.
- [x] Report each line-search trial at the largest interior residual before and after the trial.
- [x] Validate the diagnostic calculations and output without evolving a model.
- [x] Install with this checkout's `MESA_DIR` and run the package checks.

Enable in `&controls` with `RSP2_report_flux_solver = .true.`. Output goes to
the normal terminal stream with `RSP2_flux_*` prefixes and column headers.
Only one or two distinct faces are printed per trial. The surface Y boundary
row is excluded. No additional residual evaluation is performed, and the
diagnostic does not change the state, correction, residual, or tolerances.
The trace arrays are local to `star_solver:do_solver`, allocated only while
the control and RSP2 are active, and released when that solver call returns.
There is no new photo/model state and no restart-format change.

The identifiers are model, solver call, Newton iteration, line-search trial,
face k, and retry count. Each successful trial has:

- `RSP2_flux_residual`: coefficient, dt, old and new residuals, the full-row
  linear prediction, raw Newton linear residual, old/new luminosity scales,
  old/new unnormalized luminosity defects, and the current maximum tolerance.
- `RSP2_flux_Y`: Y_start, old/new Y, its column scale, raw Newton delta_Y,
  requested trial delta_Y, realized delta_Y, old dR/dY, old/new gradL, and gradT.
- `RSP2_flux_before` / `RSP2_flux_after`: L, Lr, Lc, Lt, adjacent cell pressures
  and temperatures, and their dq values. Luminosities are in erg/s, pressure
  in dyn/cm^2, and temperature in K.
- `RSP2_flux_dR`: each variable's contribution to the predicted residual
  change, summed over the three-zone stencil, in the printed variable order.

`RSP2_flux_iteration` records the trial coefficient selected by the line
search and whether the iteration passes the solver tolerances. Individual
trial records are not claims that those trials or timesteps were accepted.
`RSP2_flux_trial_failed` records evaluation errors without reporting partially
updated values as valid residuals. The call header records dynamic gradL,
u-grid flag and the matrix solver.

With the assembled, column-scaled matrix row saved before factorization,
the diagnostic computes

```math
\Delta R_k^{\rm linear}
=\sum_{j=k-1}^{k+1}\sum_i
  J^{\rm scaled}_{ki,j}
  \frac{\mathrm{solver\_dx}_{i,j}^{\rm trial}
       -\mathrm{solver\_dx}_{i,j}^{\rm before}}
       {\mathrm{x\_scale}_{i,j}},
\qquad
R_k^{\rm predicted}
=\frac{L_{{\rm scale},k}^{\rm before}}{L_{{\rm scale},k}^{\rm trial}}
 (R_k^{\rm before}+\Delta R_k^{\rm linear}).
```

The scale ratio handles the initialization evaluation at solver_iter zero,
which uses current L rather than L_start. The inward neighbor is omitted at
k=nz. `R_Newton` uses the unmodified, full Newton direction before correction
limiting; a large value identifies an inaccurate linear solve. `R_predicted`
uses the actually stored trial increment, including clipping and damping.
The difference between it and `R_after` includes nonlinear remainder, missing
or inaccurate derivatives, and errors in state evaluation; it is not by
itself a proof of a Jacobian bug. The explicit requested and realized delta_Y
and the pressure/temperature records help separate those possibilities.

The default-false toggle is also present beside the RSP2 controls in
`star/dev_cases_TDC_Pulsation/dev_TDC_RSP2_Cepheid/inlist_pulses`.
No other case settings were changed. Rebuild the case executable against this
MESA installation before running it with the new control. This task did not
build or evolve that stellar case.

Validation: the three diagnostic routine bodies were extracted verbatim into
an isolated Fortran harness and compiled with runtime bounds checks. Tests
covered full/damped/backtracked steps, unequal variable scales, neighboring
columns, the inner boundary, the initialization-to-L_start scale change,
duplicate-face suppression, a deliberately introduced nonlinear mismatch,
and a failed equation evaluation. Printed predictions, per-variable sums,
Y increments and scales were checked independently. Reporting and repeated
trials left the input solver state and the saved baseline unchanged.
The fixture and output are retained in `output/build/rsp2_flux_trace_check.f90`
and `rsp2_flux_trace_check.log`, with the source-extraction script alongside.
Control wiring/default checks, all eight independent-Y groups, and all four
remesh groups passed; `git diff --check` passed. `./install` completed with
`MESA_DIR=/Users/owner/Documents/Software/dev/test/mesa` and the package checks
passed (`output/build/rsp2_flux_diagnostics_install.log`). A final incremental
install after the boundary-copy and header polish is recorded separately in
`output/build/rsp2_flux_diagnostics_install_final.log`.

### Recorded run: a persistent inner-face residual plateau

Input: `star/dev_cases_TDC_Pulsation/dev_TDC_RSP2_Cepheid/output.txt`, restarted
from photo x00041000 at 02:04:04. The snapshot and machine-readable analysis
are in `output/review/rsp2_flux_20260919/`. Reproduce the summary with
`python3 notes/analyze_rsp2_flux_solver.py OUTPUT_TXT`.

All 840 call headers (models 41001--41840) report dynamic gradL on,
v_flag, and the banded solver. There are 839 complete tolerance
decisions before the output ends partway through model 41840; 729 passed on
iteration 10. The others passed on iterations 5--9. There are no failed
equation trials or new retries in this interval. History contains repeated
model numbers from restarts; retaining the last record per model confirms
840 accepted steps, a constant cumulative retry count of 10, and elapsed
stellar time about 0.716 day. The dt range is 35.96--101.09 seconds. The
trace ends without a termination message. This run advances; its problem is
wasted late iterations, not a demonstrated timestep collapse.

At iteration 4 the median maximum flux residual is `2.75d-8`. Its distribution
barely improves through iterations 5--10. The worst faces are overwhelmingly
347--350, not the turbulent ionization region. At iteration 10 the requested
maximum tolerance changes from `1d-8` to `1d-5`; the norm tolerance changes
from `1d-10` to `1d-8`. Most solves therefore pass because the requirements
relax, not because the flux residual suddenly drops.

For the 8575 fully recorded, coefficient-one inner-face trials at iteration
4 or later:

| Quantity | Median | Maximum |
|---|---:|---:|
| Recomputed absolute flux residual | `2.12d-8` | `1.02d-7` |
| Full-row linear prediction | `1.94d-14` | `1.94d-13` |
| Error in applying delta_Y, converted to residual units | order `1d-17` | `4.85d-17` |

The raw Newton linear residual is smaller still. Lc and Lt are exactly zero
in every one of these recorded inner-face trials. Thus neither face interpolation of w nor the negative-w repair generates
this particular late-iteration defect. These are selected worst faces, not
a census of every row of the stellar Jacobian.

`turb_support:get_TDC_dynamical_gradL` directly evaluates

```math
\nabla_{L,\mathrm{dynamic},k}
=\nabla_{L,k}\,
\frac{P_{k-1}-P_k}{-Gm_k\Delta m_{\mathrm{face},k}/(4\pi r_k^4)}.
```

At the deepest face, adjacent pressures are about `4.58d12 dyn/cm^2`, but
their difference is only `1.59d6 dyn/cm^2`, or `3.48d-7` in relative units.
The relative temperature difference is about `1d-7`. Pressure evaluation
errors or EOS derivative inconsistencies are therefore strongly amplified.
For example, model 41001, iteration 5, face 350 has:

- Requested delta_Y `-9.8327983421d-10`, applied `-9.8327983256d-10`.
- Raw Newton linear residual `8.47d-22`; stored-trial prediction `2.90d-17`.
- Recomputed flux residual `5.277326d-8`.
- Change of the actual pressure difference `-0.08203125 dyn/cm^2`.
- Actual gradL change contributes `7.0111d-8` to the flux residual;
  the assembled density-plus-temperature response contributes `1.7391d-8`.

Their difference accounts for almost all of the recomputation mismatch.
An error of only about `0.062 dyn/cm^2` in the pressure difference would
produce this residual: about `1.4d-14` of either cell pressure. That is
consistent with a few rounding units in the logarithmic EOS pressure
evaluation, but the trace alone does not prove roundoff rather than an
EOS derivative error. `micro:do_eos_for_cell/store_eos_for_cell` exponentiates
the EOS lnPgas, adds radiation pressure, and stores the total Peos.
`wrap_Peos_*` correctly uses `chiRho_for_partials` and `chiT_for_partials`.
The remaining uncertainty is the small-increment accuracy of these EOS
values relative to their derivatives, not a demonstrated missing Y column.

Subtraction of two nearby positive double-precision pressures is itself
exact under the usual factor-of-two condition; the significant uncertainty
is already in its operands. Replacing the subtraction with a logarithmic
identity or using a compensated final flux sum cannot recover those lost
digits. A stable difference evaluation would have to retain the small
thermodynamic increments through state construction and EOS evaluation, or
otherwise quantify the attainable residual accuracy. No pressure formula,
Y scale, luminosity normalization, mesh, or tolerance was changed in this audit.
A matched dynamic-gradL-off diagnostic trace would isolate the remaining
temperature/structure contribution before attempting that broader change.

### Follow-up: can the dynamic-gradL pressure evaluation be repaired?

The shared `turb_support:get_TDC_dynamical_gradL` is unchanged from committed
HEAD. Git history places the direct pressure-difference expression in
`4f4621a3c` (2026-08-31). The independent-Y RSP2 work reuses this pre-existing
TDC function; it did not introduce that subtraction. This establishes local
history only, not a newly fetched upstream comparison.

The physical expression is correct. Its numerical conditioning in thin cells
is the problem identified by this trace. The EOS path also converts the
rounded absolute lnT/lnd to base-10 logs before `eosDT_get`, then returns
absolute lnPgas. Any error introduced at those stages is amplified by the
small spatial pressure jump. A scalar 70-digit-reference calculation at a
comparable pressure/jump illustrates the limitation: direct pressure
subtraction gives a relative difference error `6.93d-10`; expm1 applied to
the already rounded logs gives `1.07d-9`; expm1 applied to the retained small
log-pressure increment gives `1.09d-15`. This is an arithmetic illustration,
not a test of the MESA EOS or a prediction of its attainable accuracy.
Values are retained in `pressure_arithmetic.json` beside the trace summary.

A repair should evaluate small EOS pressure changes from retained density
and temperature increments, keeping the residual and its endpoint derivatives
consistent. For fixed composition the exact differential is

```math
dP=P\chi_{\rho,\mathrm{partials}}\,d\ln\rho
   +P\chi_{T,\mathrm{partials}}\,d\ln T.
```

Integrating this between nearby thermodynamic states avoids differencing
large absolute pressures. A mean-derivative formula is only a quadrature
approximation and requires an error check; composition changes add their
own derivatives. Existing `chiRho_start`/`chiT_start` contain thermodynamic
EOS outputs, not necessarily the pressure derivatives used for Newton
partials, so substituting those into a pressure-increment formula without
checking consistency would be incorrect. Freezing gradL within Newton or
replacing the measured pressure gradient by hydrostatic balance changes
the selected implicit equations and is not this precision repair.

The first targeted runtime check remains the EOS pressure change versus
its derivative prediction for tiny retained increments at an affected
face. No speculative pressure approximation has been added. The evidence
strongly supports amplified evaluation noise but does not yet exclude an
EOS derivative inconsistency. A validated repair must address the origin
of that mismatch, not merely change the displayed flux residual scale.

### TDC comparison: the radiative safeguard changes the recommended repair

The user's observation that TDC behaves better prompted inspection of
`turb/public/turb.f90:set_TDC`, beyond the shared pressure-gradient helper.
There is already a specific radiative-limit safeguard at lines 210--218:

```fortran
if (conv_vel == 0d0) then
   gradT = L/L0
   Y_face = gradT - gradL
else
   gradT = Y_face + gradL
end if
```

It is present in committed `4f4621a3c` (2026-08-31). In addition, TDC solves
the face luminosity relation locally before returning the closure to the
stellar solver (`turb/private/tdc.f90:get_TDC_solution`). In the radiative
branch it can obtain Y algebraically. It does not leave an independent
global Y flux row to chase every change in the EOS-derived dynamic gradL.
In `hydro_temperature:do1_alt_dlnT_dm_eqn`, ordinary nonconvective TDC also
uses `Lrad_ad = L_ad`, whereas RSP2 takes its stored `Lr_ad`.

Consequently, the earlier explanation based only on a shared precision
problem was incomplete. The observed pressure-difference noise is real in
the RSP2 trace, but the different closure evaluation explains why TDC can
avoid this convergence symptom. The independent-Y rewrite failed to retain
TDC's numerical separation of radiative transport from the neutral gradient.
This is not a new sign/factor error in the dynamic-gradL formula.

The cleaner candidate is to use **gradT at the face as the independent
unknown**, retaining the same number of unknowns and deriving the physical
superadiabaticity:

```math
Y_{\mathrm{face},k}=\mathrm{gradT}_k-\mathrm{gradL}_k,
```

```math
R_{\mathrm{flux},k}=
\frac{\mathrm{Lrad\_coeff}_k\,\mathrm{gradT}_k
      +Lc_k(Y_{\mathrm{face}},w)+Lt_k(w)-L_k}{L_{\mathrm{scale},k}}=0.
```

The selected temperature-gradient equation consumes gradT directly. The
cell turbulent-energy equation and face/cell centering remain as before.
PII, Lc, and Source would consume derived Y through the direct algebraic
correlation. This is a coordinate change, not another unknown or a local
iteration. At fixed other state, `dR_flux/dgradT` is the same sum of positive
radiative and convective responses as the former `dR_flux/dY`.

In a radiative face where Lc and Lt vanish identically, the flux residual
then contains no dynamic gradL. Any rounding in the latter affects only
derived Y, not radiative luminosity. This is the same physical decoupling
that the TDC safeguard achieves. At an exact-zero w state that is free to
become turbulent, derivatives with respect to w must still be retained; no
new branch should freeze w or suppress emerging convection.

For an idealized perturbation of gradL with fixed thermodynamic radiative
coefficient, the current Y formulation changes the radiative flux residual
by `(Lrad_coeff/L_scale)*delta_gradL`. With gradT as the unknown that change
is exactly zero. In active convection the remaining dependence is through
`delta_Y=-delta_gradL`, which belongs to the buoyancy/enthalpy closure and
must not be discarded. This identity establishes removal of the identified
radiative noise channel, not a claim of full stellar convergence.

This candidate is preferable to broadening EOS precision work as the first
repair. It needs a consistent variable/AD definition, correction scale,
start-state storage, photo/model initialization, remesh interpolation, and
LNA transformation (`delta_Y=delta_gradT-delta_gradL`). It must preserve
current independent-Y photo/model usability through explicit conversion,
not relabel stored Y values as gradT. Simply overriding `gradT=L/Lrad_coeff`
while leaving the existing independent-Y flux row would make that row an
identity in radiative regions and leave Y unconstrained.

No such coordinate change has been made yet: it changes the user's original
choice of Y as the solver variable. The current manuscript still documents
the implemented Y formulation, not
this proposed gradT formulation. The existing observational evidence and
the EOS increment discussion above remain useful, but they do not justify
changing the shared TDC pressure-gradient helper before testing this more
direct RSP2 formulation repair.

### RSP2 alpha_t source audit, 2026-09-19

The user reports instability after enabling RSP2_alfat and requests inspection
of RSP2 source, without running MESA. This audit concerns turbulent energy
transport Lt. No Fortran or inlist
changes were made. The checks below are algebra checks, not stellar runs.

There is no automatic ramp of RSP2_alfat. `ctrls_io` assigns the control
directly, and `hydro_rsp2:compute_Lt` uses it directly. The similarly named
RSP envelope relaxation control belongs to original RSP, not RSP2.
`hydro_rsp2:set_etrb_start_vars` recomputes Lt_start with the current coefficient
before the step; it does not gradually introduce the new coefficient.

The implemented face flux is (`star/private/hydro_rsp2.f90:1144`):

```math
Lt_k = \frac{\alpha_t(4\pi r_k^2)^2\overline{\rho^2}_k\Lambda_k}
                    {\Delta m_{{\rm face},k}}
          w_{{\rm face},k}(w_k^2-w_{k-1}^2),
\qquad
w_{{\rm face},k}={\tt alfa}\,w_k+{\tt beta}\,w_{k-1}.
```

Here k=1 is the surface, w is cell centered, and Lt is face centered. Lambda
comes from the shared face mixing-length function. Current-state derivatives
are retained through automatic differentiation. For nonnegative coefficients
and w, the flux goes down the turbulent-energy gradient. The same face flux
appears with opposite signs in adjacent active turbulent-energy rows and is
included in total L by `compute_L_terms`. Its residual is
(`hydro_rsp2.f90:248`, with the time weighting at line 379):

```math
w_k^2-w_{k,\mathrm{start}}^2
+\overline{P_{t,k}}\left(\rho_k^{-1}-\rho_{k,\mathrm{start}}^{-1}\right)
+\frac{\Delta t}{\Delta m_k}
  \left(\overline{Lt_k}-\overline{Lt_{k+1}}\right)
-\Delta t\left(\mathrm{Source}_k-D_k-Dr_k+Eq_k\right)=0,
\qquad
\overline{Lt_k}=\theta_L Lt_k+(1-\theta_L)Lt_{k,\mathrm{start}}.
```

This audit does not find an explicit-current-flux update or a reversed flux
sign. However, the nonlinear solve is not robust at all small-w interfaces.
For equal cell masses, freeze the positive radius, density, mixing-length
and mass factors. The dependence of an incoming flux on the two cell
velocities is proportional to

```math
Lt=\frac{w_{\rm receiving}+w_{\rm donor}}{2}
      (w_{\rm donor}^2-w_{\rm receiving}^2),
\qquad
\frac{\partial Lt}{\partial w_{\rm receiving}}
=\frac{(w_{\rm donor}+w_{\rm receiving})
        (w_{\rm donor}-3w_{\rm receiving})}{2}.
```

The omitted positive dimensional factor does not change this derivative's
sign. Incoming flux initially increases with receiving-cell w when that w
is less than one third of the donor value. At zero w the storage derivative
2w vanishes, while this incoming-flux derivative remains positive. This is
a specific consequence of the current face-w interpolation and w-squared
energy unknown, not an assertion that the flux itself points uphill.

A two-cell check isolates precisely this implemented transport term, with
closed exterior boundaries and local source, work and damping omitted.
Use dimensionless units with the positive flux prefactor and cell masses
equal to one, initial (w_receiving,w_donor)=(0,1), dt=0.1 and theta_L=0.5.
At the initial state the residual and Jacobian are

```math
R=\begin{pmatrix}-0.05\\0.05\end{pmatrix},\qquad
J=\begin{pmatrix}-0.025&-0.075\\0.025&2.075\end{pmatrix},\qquad
J\,\delta w=-R\quad\Longrightarrow\quad
\delta w=\begin{pmatrix}-2\\0\end{pmatrix}.
```

`solver_support:Bdomain` calls `clip_so_non_negative` for RSP2 w
(`star/private/solver_support.f90:643`, helper at line 718). This correction
is therefore clipped to (0,0), leaving the residual unchanged. Nevertheless,
the same centered equations have the positive solution
(0.22785876226600094, 0.9736941945285522), with residuals below 2e-17 in a
Python algebra check. Thus a valid positive root can exist while Newton and
the existing clipping fail to reach it. This counterexample is not a
reproduction of the fully coupled stellar calculation; local source terms
and the other stellar equations can change its Newton direction.

The startup path does not cover this case:

- `RSP2_adjust_vars_before_call_solver` skips every cell with nonzero w_start
  and every cell without positive local buoyant driving
  (`hydro_rsp2.f90:1223`, `hydro_rsp2.f90:1230`). It does not include Lt inflow
  in its predicted w. `struct_burn_mix:set_xh` calls it and then packs the
  updated w into the solver state, so this is an active code path.
- The positive-branch residual is explicitly limited to RSP2_alfat=0
  (`hydro_rsp2.f90:284`). That restriction cannot simply be removed: an
  incoming turbulent flux does not have local w as a common factor.
- The same face-w dependence and local-only startup restriction are present
  in the checked-out HEAD version, before the current uncommitted rewrite.
  This comparison does not establish their history on a remote branch.

These checks establish a missing treatment of turbulence entering a quiet
cell. They do not identify the cells or terms responsible for the user's
reported crash. A repair must address transport-driven onset and bound
handling while preserving the full coupled residual, shared fluxes and
centered time integration. A transport-aware positive initial guess is a
candidate for further derivation, not an implemented or validated repair.
No change to L_theta, an arbitrary turbulence floor, or a coefficient ramp
has been introduced or established as a remedy.

### RSP2 transport startup predictor implemented, 2026-09-19

On the user's instruction, extended
`hydro_rsp2:RSP2_adjust_vars_before_call_solver` to supply a positive initial
w guess for cells whose incoming turbulent energy over the step exceeds
their starting turbulent energy. This includes both exactly zero and tiny
positive w. It is an initial-guess change only. The turbulent-energy residual,
physical Lt and Lc, Jacobian assembly, luminosity time weighting, negative-w
clipping, and accepted starting state are unchanged. No MESA compilation,
installation, or stellar run was performed.

The selection uses the saved, signed flux divergence:

```math
{\tt detrb}_k=\frac{\Delta t}{\Delta m_k}
 (Lt_{k+1,\mathrm{start}}-Lt_{k,\mathrm{start}}),
\qquad
{\tt detrb}_k>w_{k,\mathrm{start}}^2.
```

This selects weak turbulence relative to the incoming energy; it introduces
no absolute w threshold. Net loss, zero transport, zero timestep, and
RSP2_alfat=0 do not activate the new predictor. Forced nonturbulent cells
are excluded using the turbulent-energy equation's existing mask. At the
innermost cell the missing inner face flux and its derivative are zero.

The first candidate, sqrt(w_start^2 + detrb), was rejected by an algebra
check before completion: for an outer receiving cell of mass 1e-4 next to
a donor of mass 1, the mass-weighted face velocity initially gives very
little weight to the donor. Increasing receiving-cell w rapidly increases
that face velocity. A flux-only estimate still lands on the wrong Newton
branch. The implemented predictor therefore uses the existing AD derivative
of Lt with respect to the receiving cell's w as well as the flux itself.

Before bounding the guess, it solves the following quadratic exactly, with
neighboring states and the other current variables held fixed:

```math
w_{k,\mathrm{guess}}^2=w_{k,\mathrm{start}}^2
+\frac{\Delta t}{\Delta m_k}\left[
 (1-\theta_L)(Lt_{k+1,\mathrm{start}}-Lt_{k,\mathrm{start}})
+\theta_L\left\{Lt_{k+1}-Lt_k+
 \left(\frac{\partial Lt_{k+1}}{\partial w_k}
       -\frac{\partial Lt_k}{\partial w_k}\right)
 (w_{k,\mathrm{guess}}-w_k)\right\}\right].
```

The storage term remains quadratic, while only the transport term is
linearized for this initial estimate. Theta_L is read using exactly the
same controls as the physical turbulent-energy row. In the helper,
`detrb_dw` is the time-weighted derivative of the incoming specific energy
with respect to w, and `etrb` is the remaining constant term in the quadratic
`w_guess^2 - detrb_dw*w_guess - etrb = 0`. The positive root is evaluated
without subtractive cancellation. Invalid coefficients leave the current
guess alone. The new guess is bounded above by the largest neighboring
starting w and only replaces a smaller current guess. This bounds the initial
estimate, not Lt or the converged w; subsequent Newton iterations still use
the original nonlinear residual without this bound.

Source references and integration:

- `star/private/hydro_rsp2.f90:1196`: save the already computed Lt AD value
  in the existing `Lt_ad` cache when saving Lt_start. This keeps the value
  and derivatives consistent for the first prediction, including a restart
  where RSP2_alfat has changed.
- `star/private/hydro_rsp2.f90:1218`: transport predictor before the existing
  local buoyancy predictor. The outer face uses derivative `i_w_00`; the
  inner face uses `i_w_m1`, because cell k is the outer neighbor of face k+1.
  The loop reads fixed flux caches and w_start values and writes only its
  own cell's w, so it does not depend on OpenMP loop order.
- The original local-buoyancy loop body and its effective bounds are
  unchanged. Transport prediction also applies with nonzero source_seed;
  that control continues to bypass only the original local predictor.
- `struct_burn_mix:save_start_values` runs before `set_xh`; `set_xh` packs
  the predicted w, and `do_solver` records its difference from xh_start in
  solver_dx. There is no write to w_start, Lt_start, or xh_start in the helper.
  A retry saves the restored starting state and evaluates the predictor
  using the retry timestep.
- Existing `RSP2_report_adjust_w` also reports these updates with the label
  `RSP2_adjust_vars_before_call_solver Lt w`. No new control was introduced.

Validation in `notes/check_rsp2_transport_startup.py`:

- The original equal-mass, zero-w counterexample stalls without prediction
  and converges in three Newton updates with it, to the same positive root
  recorded above.
- The unequal-mass example that defeated the flux-only estimate also
  converges with the derivative-aware estimate.
- All 672 isolated two-cell transport cases converged to the original
  nonlinear residual tolerance, with a maximum of five Newton updates.
  Cases cover mass ratios 1e-4 through 1e4, both flux directions, equal and
  mass-weighted face interpolation, zero/1e-14/1e-8 initial receiving w,
  theta_L=0.5 and 1, and four timestep scales. The converged roots conserve
  the two-cell turbulent energy to the checked floating-point tolerance.
- Independent finite differences checked the coupled transport Jacobian;
  checks also cover unchanged starting arrays, zero-step and uniform-state
  behavior, and the neighboring-value bound for a long-step guess.
- Source comparison against the pre-edit snapshot confirms that all physical
  fluxes and residuals, and the original local-startup loop body, are byte
  unchanged. Boundary/mask checks include a one-cell mesh and an active
  last cell with no inner stored face. `git diff --check` passes.

The script checks reduced transport equations, not MESA's full coupled
Jacobian or source/work terms. This implements and checks the identified
startup correction; it does not establish that the reported stellar crash
is cured. Full stellar convergence and pulsation behavior remain untested.
The pre-edit source snapshot is in
`output/review/rsp2_alfat_startup_20260919/hydro_rsp2.before.f90`.

### Transport startup predictor installed, 2026-09-19

The user authorized installation. `./install` returned exit code 0 with
`MESA installation was successful`; it compiled the changed hydro_rsp2 and
passed the installer package checks. Both installation and the subsequent
Cepheid executable relink explicitly exported
`MESA_DIR=/Users/owner/Documents/Software/dev/test/mesa` after initializing
`MESASDK_ROOT=/Applications/mesasdk`, with GYRE_DIR unset and NPROCS=10.
The installer records this MESA_DIR in build.log.

Ran `make -j10 -W build/src/run.o all` in
`star/dev_cases_TDC_Pulsation/dev_TDC_RSP2_Cepheid` to force the executable
relink against this checkout's updated libraries; it returned exit code 0.
The new transport-predictor diagnostic string was verified in hydro_rsp2.o,
libstar.a, and the case's build/bin/star. The object, library, and executable
timestamps are newer than the source. Case inlists were not edited and no
Cepheid evolution or standalone diagnostic executable was run.

Logs and source/binary hashes:

- `output/build/rsp2_transport_startup_install_20260919.log`
- `output/build/rsp2_transport_startup_case_build_20260919.log`
- `output/build/rsp2_transport_startup_install_20260919.json`

This supersedes the uncompiled/uninstalled status above. Full stellar
convergence and pulsation behavior still require the user's model test.

### work2 alpha_m and conservative-work Jacobian audit, 2026-09-19

The user reports that `star/work2` cannot get through startup with
RSP2_alfam=0.25 and asks whether omitted derivatives at shift_p1 are
responsible. Source inspection confirms incomplete viscous Jacobians in
the v_flag path. This is distinct from the alpha_t startup predictor above.
No source changes, installation, or MESA runs were made during this audit.

The supplied terminal paste contains four starts of the 15-Msun pre-MS
case. Each completes the pre-MS preparation through relaxation model 34,
then enables RSP2 and prints `v_flag T`. No accepted RSP2 step follows in
the supplied intervals. The last attempt has solver calls 336--343 after
the switch, seven retries, and a final retry timestep log10(dt/yr) of
-7.107209969647868. Several attempts reach 250 iterations. The persistent
late maximum residual is dlnE_dt in zone 2; one 250-iteration attempt ends
at 1.6232864e-3 against a 1e-4 maximum-residual tolerance. This is first-step
stagnation, not evidence of hundreds of accepted RSP2 models.

The current inlist sets energy_eqn_option='dedt' and does not override
use_P_d_1_div_rho_form_of_work, whose default is false. Although the current
inlist says new_v_flag=false, `set_flags:set_RSP2_flag` enables v_flag when
u_flag is absent (`star/private/set_flags.f90:393`). The paste confirms that
this happened. The current file has subsequently been set to RSP2_alfam=0;
the paste does not print alpha_m for each attempt, so the user's report
supplies the identification of the failing alpha_m=0.25 case.

#### Missing inner second-neighbor derivative in total energy

`hydro_energy:setup_sources_and_others` includes the following viscous source
for RSP2 v_flag in total-energy form (mass corrections omitted here for
clarity):

```math
Eq_k+\frac12\left(\overline v_k Uq_k+
                  \overline v_{k+1}Uq_{k+1}\right),
\qquad \overline v_j=\frac{v_j+v_{j,\mathrm{start}}}{2},
\qquad
Uq_{k+1}=\frac{4\pi(\Chi_k-\Chi_{k+1})}
 {r_{k+1,\mathrm{work}}\Delta m_{\mathrm{face},k+1}}.
```

But Chi_{k+1} contains the shear
`v(k+1)/r(k+1) - v(k+2)/r(k+2)` as well as its mixing length and radius
factors. Thus the energy row has real k+2 dependencies. At
`star/private/hydro_energy.f90:300`, the inner acceleration is formed as

```fortran
Uq_p1 = shift_p1(compute_Uq_face(s, k+1, ierr))
```

`auto_diff_support:shift_p1` (`star/private/auto_diff_support.f90:35`) retains
the value, maps old m1 to 00 and old 00 to p1, and discards old p1. That
discarded slot belongs to physical zone k+2 here. The matrix receives no
derivative for it: the hydro row storage supports only k-1, k and k+1.
This is an incomplete Jacobian, not an omitted residual-value work term.
The same expression is already present in the checked-out HEAD version of
hydro_energy; it was not introduced by the latest alpha_t predictor.

For simple P d(1/rho) or eps_grav, include_dke_dt is false and the explicit
v*Uq work term is absent (`hydro_energy.f90:94`). This particular energy-row
coupling is consequently absent. That does not establish that all remaining
viscosity derivatives fit the current stencil.

#### Additional outer dependency from the shared cell mixing length

The current RSP2 cell stress uses `get_TDC_mixing_length_cell`, whose actual
definition is (`star/private/tdc_hydro.f90:161`)

```math
\Lambda_{\mathrm{cell},j}
 =\tfrac12(\Lambda_{\mathrm{face},j}
                 +\Lambda_{\mathrm{face},j+1}).
```

The outer face scale height depends on the EOS states of cells j-1 and j
(`star/private/star_utils.f90:4087`). Consequently Chi_{k-1} depends on
thermodynamic variables in cell k-2 once shear is nonzero. However,
`compute_Uq_face` shifts Chi_{k-1} with shift_m1
(`star/private/hydro_rsp2.f90:779`), dropping those outer derivatives.
This affects the momentum Jacobian even with simple PdV work.

There is also a nested-shift loss inside the nominal energy-row stencil:
compute_Uq_face(k+1) first shifts Chi_k into the k+1 frame and drops Chi_k's
k-1 dependency. The outer shift_p1 back to the energy row cannot recover it,
although k-1 is a representable column of that energy row.

This additional thermodynamic dependence follows from our uncommitted switch
to the shared face-averaged cell length. The HEAD RSP2 stress used a local
P/(rho*g_cell) cell scale height instead. Retaining the requested shared
scale-height functions requires handling the resulting derivatives, not
silently treating those functions as cell-local.

For the u_flag RSP2 formulas with the default two-cell face states, each
face stress depends on its two adjacent cells. Shifting the inner stress
there fits the three-cell stencil for both Eq_cell and Uq_dm_cell. The
specific v_flag omissions above therefore do not apply to that arrangement.
This observation is not a blanket Jacobian certification for every face
reconstruction option.

An independent algebra check evaluated the actual stress, force and energy
expressions on seven cells using an ideal-gas P=rho*T state, equal dm, HSE
face heights, alpha_m=0.25 and mixing_length_alpha=2. Full derivatives from
complex steps agree with centered finite differences. Emulating the native
three-cell shift operations instead gives the following representative
unnormalized derivatives for row k=4:

| Derivative | Full expression | After the implemented shifts |
| --- | ---: | ---: |
| viscous energy source with respect to v(k+2) | 126.881258388 | 0 |
| Uq(k) with respect to lnT(k-2) | 454.274126571 | 0 |
| viscous energy source with respect to lnT(k-1) | 125.394262054 | 179.989881608 |

These numbers demonstrate the stencil loss; they are not measurements from
the stellar run. Previous telescoping-conservation checks tested residual
values, not the complete spatial Jacobian, and did not detect these defects.

The conservative equation itself remains a valid target. A complete repair
must retain or exactly eliminate the additional couplings while preserving
the momentum/energy transfer and shared scale heights. Replacing the work
form or dropping the v*Uq term is not such a repair. The source defects are
confirmed; the supplied trace alone does not establish how much of this
particular startup failure each defect causes. No repair is implemented yet.

### TDC comparison and scope of a cell-length rewrite, 2026-09-19

The user asks why the cell-length dependence would be a problem when TDC
works. Direct comparison confirms that this omission is also present in
the current TDC v_flag eddy-viscosity path when TDC_alpha_M is nonzero:
`tdc_hydro:compute_Chi_div_w_cell` calls the same face-averaged
get_TDC_mixing_length_cell (line 201), and `compute_tdc_Uq_face` shifts the
outer cell stress with shift_m1 (line 479). RSP2 inherited this dependence
by adopting the shared cell helper. TDC is not evidence that the resulting
Jacobian is complete. The code and observed convergence must be distinguished.
These thermodynamic stress derivatives are proportional to the background
shear, so they vanish at a static zero-shear state and can be weak in other
regimes; an approximate Jacobian can still converge.

The explicit dependency is
`Uq(k) -> Chi(k-1) -> Lambda_face(k-1) -> EOS(k-2)`.
It is real dependence of the current discrete formula, not an error in
relabeling a derivative that could simply be assigned to another slot.
Expanding the same expression algebraically cannot remove it. The inner
total-energy work coupling discussed above is a separate issue.

A compact way to remove this particular outer thermodynamic dependence is
to define the shared cell mixing length from the cell's own EOS state,
cell gravity and cell radius, applying the same scale-height options and
get_mlt_mixing_length law as appropriate. Both TDC and RSP2 could then use
that shared cell definition; face lengths would retain their face definition.
This is a deliberate change in cell discretization from averaging two face
lengths, not an algebraically identical rewrite. It should not be silently
presented as preservation of the existing face-averaged formula.

Preserving that exact face-averaged formula instead requires retaining the
wider derivative stencil or introducing an algebraic stress unknown with
its own constitutive residual. The latter can keep nearest-neighbor blocks
while preserving the eliminated equations, but requires solver/state plumbing
and is not a small helper-only edit. Neither option has been implemented in
this clarification. A local cell-length change alone would not repair the
separate k+2 velocity dependence of the conservative energy work term.

### Shared cell-local mixing length and static LNA update, 2026-09-19

Implemented the cell-length correction described above. This supersedes the
face-average definition and the "not implemented" status in the preceding
audit. Only `tdc_hydro.f90` and `star_LNA_turbulence_closures.f90` required
new source edits. No MESA compilation, installation, or model run was made
for this change.

#### Cell definition, options and boundaries

`tdc_hydro:get_TDC_mixing_length_cell` (lines 161--188) now evaluates

```math
g_{\mathrm{cell},k}=\frac{g_k+g_{k+1}}2,\qquad
g_j=\frac{\mathrm{cgrav}_j m_{\mathrm{grav},j}}{r_j^2},\qquad
H_{P,\mathrm{cell},k}^{\mathrm{HSE}}
 =\frac{P_{\mathrm{eos},k}}{\rho_k g_{\mathrm{cell},k}}.
```

P and rho use `wrap_Peos_00` and `wrap_d_00`, retaining the cell's EOS
partials. The gravity average uses `wrap_r_00` and `wrap_r_p1`. At the full
centre, the inner gravity is zero, without division by zero. At an excised
boundary, it is cgrav(nz)*M_center/R_center**2; R_center and M_center remain
fixed boundary values. The outer face uses m_grav, including the existing
mass corrections. These heights and lengths have units of cm.

With harmonic_dissipation_length_beta <= 0 and alt_scale_height_flag true,
the selected height is

```math
H_{P,\mathrm{cell},k}
 =\min\left(H_{P,\mathrm{cell},k}^{\mathrm{HSE}},
            \frac{\sqrt{P_{\mathrm{eos},k}/\mathrm{cgrav}_k}}{\rho_k}\right).
```

Otherwise it is the HSE height. The positive harmonic option continues to
take precedence over the alternate height, matching the face length policy.
The existing `star_utils:get_mlt_mixing_length` then evaluates

```math
\Lambda_{0,k}=\alpha_{\mathrm{MLT}}H_{P,\mathrm{cell},k},\qquad
r_{\mathrm{cell},k}=\left[\frac{r_k^3+r_{k+1}^3}{2}\right]^{1/3},
\qquad
\Lambda_{\mathrm{cell},k}=
\begin{cases}
 \Lambda_{0,k}, & \beta_h\le0,\\
 [\Lambda_{0,k}^{-1}+(\beta_h r_{\mathrm{cell},k})^{-1}]^{-1},&\beta_h>0.
\end{cases}
```

The AD cell radius matches the volume midpoint in `star_utils:set_rmid`.
There is no half-length special case at the innermost cell. Consequently
this is a different cell discretization, not a rearrangement of the former
face average. Changes need not be small near a centre, steep gradients, or
an alternate-height transition.

At fixed gravitating masses and cgrav, the differentiated HSE branch is

```math
\delta\ln H_{P,\mathrm{cell},k}^{\mathrm{HSE}}
 =\delta\ln P_{\mathrm{eos},k}-\delta\ln\rho_k
  +2\frac{g_k\delta\ln r_k+g_{k+1}\delta\ln r_{k+1}}{g_k+g_{k+1}}.
```

The alternate branch has delta ln H = 0.5 delta ln P - delta ln rho.
For the harmonic branch,

```math
\delta\ln\Lambda_{\mathrm{cell},k}
 =\frac{\Lambda_{\mathrm{cell},k}}{\Lambda_{0,k}}\delta\ln H_{P,\mathrm{cell},k}
 +\frac{\Lambda_{\mathrm{cell},k}}{\beta_h r_{\mathrm{cell},k}}
    \delta\ln r_{\mathrm{cell},k},\qquad
\delta\ln r_{\mathrm{cell},k}
 =\frac{r_k^3\delta\ln r_k+r_{k+1}^3\delta\ln r_{k+1}}{r_k^3+r_{k+1}^3}.
```

Terms from the fixed inner radius have zero variation. The cell length has
only EOS(k), r(k) and r(k+1) dependencies. Therefore shifting Chi(k-1) into
the momentum row no longer drops an EOS(k-2) dependency introduced by the
length. The former nested-shift loss of EOS(k-1) in the inner Uq contribution
also disappears. The conservative-work dependence on k+2 remains.

#### Consumers, reconstruction and LNA

- TDC v-grid viscosity uses the shared helper in `compute_Chi_div_w_cell`
  (`tdc_hydro.f90:213`). Its face convective velocity and its averaging into
  the cell stress remain unchanged. TDC face convection and u-grid viscosity
  retain their face length.
- RSP2 uses the shared cell length in `compute_D_div_w` (line 582),
  `compute_Dr_div_w` (line 626), `compute_Chi_div_w_cell` (line 669), and the
  existing local startup estimate (line 1245). Thus its damping coefficients
  change along with the v-grid stress. RSP2 w remains cell centred.
- Both `rsp2_chi_coefficient_for_star_LNA` (line 1207) and
  `tdc_chi_coefficient_for_star_LNA` (line 1223) now call the shared cell
  helper. The duplicate LNA-only face average was removed. RSP2 LNA damping
  already calls nonlinear `compute_D` and `compute_Dr` through its wrappers
  (lines 832 and 845), so their AD derivatives also use the new definition.
- `get_TDC_Hp_face` and `get_TDC_mixing_length_face` are byte-for-byte
  unchanged by this patch. They still use reconstructed face states when
  use_face_reconstruction is true. Cell lengths use cell EOS quantities
  in either mode; reconstructing a face does not change that cell definition.
  use_rsp_form_of_scale_height still controls ordinary face interpolation;
  at a cell there is only one P/rho value, so no such interpolation choice
  is needed. No claim is made that every reconstruction derivative has been
  independently validated by this patch.

At a static zero-shear background, the LNA viscosity linearization is

```math
\delta\Chi_k=
 \left[\frac{16\pi\alpha_m}{3\Delta m_k}\rho_k^2
       \frac{r_k^6+r_{k+1}^6}{2}\Lambda_{\mathrm{cell},k}w_{\mathrm{stress},k}\right]_0
 \left(\frac{\delta v_k}{r_k}-\frac{\delta v_{k+1}}{r_{k+1}}\right).
```

Perturbations of the coefficient multiply zero background shear. Taking
Lambda_cell%val for this static coefficient is therefore consistent; RSP2
damping, which is nonzero in equilibrium, retains the full AD length.

#### Conservation and checks

The patch changes the shared stress coefficient, without changing Eq or Uq
quadrature, time centering, masks, or the energy-equation branches. For the
paired v-grid stress and work shear, summation by parts gives

```math
\sum_k\Delta m_k Eq_k
 +\sum_j\Delta m_{\mathrm{face},j}\overline v_j Uq_j
 =-4\pi\Chi_N\frac{v_{\mathrm{inner,work}}}{r_{\mathrm{inner,work}}},
```

where the right side is inner-boundary work and is zero at a full centre or
fixed zero-velocity boundary. This identity holds for any common Chi and
does not require the old face-averaged cell length. It is not a certification
of unrelated TDC explicit-w or boundary-mask policies, which were not changed.
RSP2 turbulent/gas damping transfers retain the same equal-and-opposite
values, even though the damping length changes. The u-grid stress and work
functions are unchanged. Both dedt work forms and eps_grav keep their existing
selection logic.

`notes/check_tdc_cell_mixing_length.py` is an independent Python algebra and
source check, not execution of MESA's Fortran AD code. Results:

- 72 length cases covering HSE, alternate and harmonic branches, harmonic
  length ratios from 1e-8 to 1e8, full centre and excised/interior geometry:
  complex-step versus analytic log derivatives agree within 8.9e-16.
- 16 combinations of boundary geometry, length options and RSP2/TDC w
  placement: no stress dependency outside cells k and k+1; the static LNA
  coefficient agrees with the simultaneous perturbation of geometry, EOS,
  turbulent velocity and radial velocity at zero background shear.
- 100 unequal-mass, mass-corrected viscous-power cases with arbitrary stresses
  and inner-boundary work: the identity above closes within 1.0e-15 normalized
  error. This checks residual values, not the unresolved energy-row Jacobian.
- Source checks confirm unchanged face and Eq/Uq routines and the shared LNA
  calls. Fortitude passed both modified Fortran files; git diff --check passed.

The manuscript's cell-length formula and TDC LNA explanation were updated
and its 32-page PDF rebuilt. Pages 12--14 and 24 were rendered and visually
reviewed; the changed equations and surrounding row boxes fit without clipping.
Source snapshots and build output are in
`output/review/cell_mixing_length_20260919/`. Subsequent installation is recorded
below; convergence of the user's models remains untested for this change.

#### Installation, 2026-09-19

At the user's explicit request, ran `./install` after initializing
`MESASDK_ROOT=/Applications/mesasdk` and then explicitly setting
`MESA_DIR=/Users/owner/Documents/Software/dev/test/mesa`, with NPROCS=10,
OMP_NUM_THREADS=8, and GYRE_DIR unset. Installation and the standard package
checks passed with exit status 0. Both changed Fortran sources were compiled;
their resulting objects were checked byte-for-byte against the corresponding
members of `build/star/lib/libstar.a`.

`star/work2` and `dev_TDC_RSP2_Cepheid` were relinked successfully before the
user clarified the standing preference: **the user handles all case relinking;
future install requests must stop at the MESA installation.** Neither user
model was run.

Logs and source/object hashes are in
`output/build/shared_cell_length_install_20260919.log` and
`output/build/shared_cell_length_install_20260919.json`. The two completed
case-build logs are `shared_cell_length_work2_build_20260919.log` and
`shared_cell_length_cepheid_build_20260919.log` in the same directory.

### Signed PII closure cleanup, 2026-09-19

At the user's request, `hydro_rsp2:compute_PII_from_Y` now contains only the
signed algebraic correlation and its shared face-length calls. RSP2 LNA
already reads the same PII_ad for luminosity and calls the nonlinear source;
both consumers were verified unchanged. Source, damping, energy/work forms,
time centering, cell/face placement, and restart formats retain their existing
definitions. The independent Y and L unknowns are unchanged.

The TDC control and its implementation in `turb/private/tdc_support.f90` and
`turb/public/turb.f90` are retained byte-for-byte. Its defaults documentation
now states that it applies only to TDC. The obsolete entry was removed from
the RSP2 Cepheid inlist, and the RSP2 flux-diagnostic header no longer reports
that TDC-only setting. The three obsolete cap-related profile columns were
removed; the other 14 RSP2 extra columns, including start values, face w and
Source_div_w, remain. This changes diagnostics, not stored model/photo state.

The active notes and manuscript now give the signed PII relation and its
product-rule perturbation. Superseded experiment narratives and two analysis
scripts were removed from the active notes after snapshotting them in the
review directory cited at the top. External reference papers are unedited.
The manuscript retains its energy-form, work, grid-placement and boundary
details; the defining PII equation and its perturbation are kept on one page.

Validation completed without compiling or executing MESA:

- All eight groups in `check_rsp2_independent_Y.py` pass. The PII derivative
  check covers four unequal face weights and nine signed Y values from
  -1e8 through zero to +1e8, using complex-step derivatives. The same test
  retains source/flux staggering, work-form, conservation and remap checks.
- `check_tdc_cell_mixing_length.py` passes, including static LNA consistency.
- Fortitude passes all three newly edited Fortran files; git diff --check
  passes. Updated Python scripts parse successfully.
- Profile-column assignments match their advertised count of 14. Both TDC
  source files and both LNA consumer modules match their pre-cleanup snapshots.
- The 32-page PDF rebuild has no unresolved references or overfull boxes.
  All pages were inspected in an overview and revised pages 24--26 were
  rendered at higher resolution and reviewed. The pre-existing underfull
  paragraph on the cell-length page remains visually acceptable.

No installation, case relinking, or user model run was performed for this
cleanup. The user handles case relinking. Full runtime robustness and the
separate k+2 conservative-work Jacobian issue remain to be verified/addressed.
The additional two turbulence equations remain deferred; no solver slots
were allocated or repurposed for them.

### Three-equation literature and initializer review, 2026-09-19

At the user's request, `notes/rsp2_three_equation_local_proposal.md` records
the two extra entropy-moment equations with local relaxation in place of
their nonlocal transport, their two published closure coefficients, and
their relation to current RSP2 variables. Original Flaskamp, Braun, and
Ahlborn sources were checked. It distinguishes the old one-equation Dr
sink, actual entropy-gradient driving, and omitted mean-strain terms.
It also derives an analytic local equilibrium and a flux-preserving
snapshot seed, including their limits and a continuous realizability concern.
Standalone algebra checks passed; no MESA source or manuscript changes,
compilation, installation, relinking, or model runs were made for this review.
The follow-up staggering review recommends both new moments on faces for
the local version, retaining cell w^2. The proposal records face energy
reconstruction, cell buoyancy exchange, and boundary requirements. This
remains a design recommendation; no solver state was changed.
