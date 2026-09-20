# RSP2 Flaskamp temperature-gradient reformulation

Historical formulation: the active RSP2 implementation now solves independent
face Y and a separate luminosity balance, as documented in
[rsp2_independent_Y_implementation.md](rsp2_independent_Y_implementation.md).
The eliminated-gradient derivations below do not describe the current solver.


Status: source implementation complete on 2026-08-28. Scoped Fortitude,
`git diff --check`, and a clean `./install` pass. Partials checks and model
tests remain to be run.

## Motivation

The previous interior RSP2 equation imposed

\[
R_L = \frac{L_r+L_c+L_t-L}{L_{\rm scale}}=0.
\]

When a TDC photo is changed to RSP2, the stellar state retains its temperature,
density, luminosity, and radial velocity. RSP2 initializes its turbulent
velocity from `mlt_vc`, but the resulting RSP2 flux partition does not generally
satisfy the direct luminosity equation. In the model-6000 restart, the first
RSP2 solve had a maximum `equL` residual of 139.56 in zone 4. Later iterations
developed large `w` corrections and outer momentum residuals.

The direct flux residual was poorly conditioned because a small change in the
superadiabatic gradient can produce a large change in luminosity. The Flaskamp
one-equation formulation instead eliminates the temperature gradient from the
local luminosity balance and uses the ordinary MESA temperature-structure row.

## Algebraic closure

Let

\[
L_r=K_r\nabla_T,
\qquad
L_c=K_c\mathcal{Y},
\qquad
\mathcal{Y}=\nabla_T-\nabla_{L,*},
\]

where `K_r` is the radiative luminosity per unit gravity-normalized temperature
gradient and `K_c` is the linear RSP2 enthalpy-flux coefficient. The neutral
gradient \(\nabla_{L,*}\) is selected as described below. Including the
nonlocal turbulent luminosity gives

\[
L-L_t=K_r\nabla_T+K_c(\nabla_T-\nabla_{L,*}).
\]

The local temperature gradient is therefore

\[
\boxed{
\nabla_T=
\frac{L-L_t+K_c\nabla_{L,*}}
     {K_r+K_c}
}.
\]

The identity

\[
L=L_r+L_c+L_t
\]

then follows algebraically. It is not imposed as a separately normalized
Newton row.

The interior `equL` row now uses the existing MESA temperature equation in
`hydro_temperature::do1_dlnT_dm_eqn`,

\[
R_T=
\Delta m
\left(\frac{d\ln P}{dm}\right)_{\rm QHSE}
\nabla_T-\Delta\ln T=0.
\]

The surface equation is selected by `RSP2_use_L_eqn_at_surface`.

## Dynamical neutral gradient

RSP2 should honor `TDC_use_dynamical_gradL` so that RSP2 and TDC use the same
neutral-gradient coordinate.

When the control is false,

\[
\nabla_{L,*}=\nabla_L.
\]

When the control is true, use
`turb_support::get_TDC_dynamical_gradL` with full automatic derivatives:

\[
\nabla_{L,*}=f\nabla_L,
\qquad
f=\frac{\Delta P_{\rm actual}}{\Delta P_{\rm QHSE}}.
\]

The helper uses the instantaneous hydrodynamic state, mass corrections, and
the same rotation correction as the momentum equation. It returns the
ordinary \(\nabla_L\) at a boundary, without hydrodynamics, or when the QHSE
pressure interval is invalid.

The pressure-gradient factor is not velocity-time-centered. Its discrete
definition is

\[
f_k=
\frac{P_{k-1}-P_k}
{-g_k\overline{\Delta m}_k/(4\pi r_k^2)},
\]

with every quantity evaluated on the current grid. This is the spatial RSP
factor, not a time-discretized momentum term. Mixing current pressure and
radius with `Peos_start` and `r_start` is also invalid immediately after a
remesh because the start arrays can still describe the old grid. The TDC
envelope reconstruction enforces the current-grid QHSE pressure jump, so this
form gives \(f_k\simeq1\) immediately after that reconstruction.

With the control enabled, the RSP2 closure is

\[
\mathcal{Y}=\nabla_T-f\nabla_L,
\]

and the algebraic gradient is

\[
\boxed{
\nabla_T=
\frac{L-L_t+K_c f\nabla_L}
     {K_r+K_c}
}.
\]

This is the same coordinate used by the current TDC dynamical-`gradL` option.
The factor belongs only on the neutral-gradient intercept. It must not multiply
`gradr`, `K_r`, the supplied luminosity, or the complete radiative term. In
particular, the following forms are incorrect:

\[
L_r=fK_r\nabla_T,
\qquad
\nabla_T=f\frac{L-L_t+K_c\nabla_L}{K_r+K_c}.
\]

For a Ledoux calculation, the implementation target is exact parity with the
current TDC path: pass the complete stored
\(\nabla_L=\nabla_{\rm ad}+B\) through
`get_TDC_dynamical_gradL`. This gives `f*gradL` wherever that helper applies.

MESA smooths the pressure-normalized `brunt_B` before storing the composition
term in `gradL`. Multiplying that stored value by a sharply varying local
\(f\) is not identical to forming and smoothing the gravity-normalized
composition term. The pressure-inversion fallback in `get_brunt_B` also
already uses a QHSE pressure interval. Full Ledoux support therefore still
requires a separate, consistently smoothed \(B_{\rm QHSE}\), with

\[
\nabla_{L,*}=f\nabla_{\rm ad}+B_{\rm QHSE}.
\]

The current implementation retains `f*gradL` and documents this limitation.

### Radiative Ledoux branch

A restart from model 1000 exposed a separate cancellation at a sharply
resolved composition interface. Near face 570, the model had
\(dq\simeq7.1\times10^{-11}\), \(f\simeq27\), zero convective velocity, and
an inferred composition contribution \(B\simeq5.2\times10^4\). Neighboring
faces reached \(B\simeq2.9\times10^5\). The exact TDC solution there is

\[
v_c=0,
\qquad
\nabla_T=\frac{L}{L_0}.
\]

`set_TDC` previously recovered this gradient as

\[
\nabla_T=\left(\frac{L}{L_0}-\nabla_{L,*}\right)+\nabla_{L,*}.
\]

This subtracts and then adds the large Ledoux term in both the value and its
Jacobian. The resulting loss of precision appeared as repeated `equL`
failures near the composition interface. When the solved convective velocity
is zero, `set_TDC` now returns `L/L0` directly and then defines
`Y_face = gradT - gradL`. Convective faces retain the existing reconstruction.

## Implemented evaluation order

`hydro_rsp2::set_RSP2_vars` assembles the RSP2 state in this order:

1. Form `L_t` from the current turbulent-velocity field.
2. Form `K_r` from the current face thermodynamic and opacity state.
3. Form the linear `K_c` from the current RSP2 turbulent velocity and face
   thermodynamic state.
4. Set `gradL_star = gradL`, or call `get_TDC_dynamical_gradL` when
   `TDC_use_dynamical_gradL` is true.
5. Evaluate the algebraic `gradT_ad` expression above with full AD derivatives.
6. Use the same `gradT_ad` and `gradL_star` in `Y_face`, `PII`, the turbulent
   source, `L_c`, and `L_r`.
7. Route the interior `equL` row through
   `hydro_temperature::do1_dlnT_dm_eqn` while retaining the RSP2 surface
   luminosity row.

The generic MESA temperature row cannot be selected without first replacing
its input `s% gradT_ad` by the RSP2 algebraic closure. Otherwise the row would
use the ordinary MLT or TDC gradient and would not represent RSP2.

`compute_RSP2_gradT` is a private helper in `hydro_rsp2`. It does not add a
Newton variable, residual row, restart field, or wrapper.

## Controls and boundaries

`RSP2_use_RSP_eqn_for_Y_face` was removed. RSP2 now has one definition,

\[
\mathcal{Y}=\nabla_T-\nabla_{L,*},
\]

for the turbulent source, `PII`, and convective luminosity.

`RSP2_use_mass_interp_face_values = .true.` applies the existing mass weights
to every RSP2 interpolation from adjacent cell centers to their shared face,

\[
X_f=\frac{\Delta m_{k-1}X_k+\Delta m_kX_{k-1}}
          {\Delta m_{k-1}+\Delta m_k}.
\]

`hydro_rsp2::get_RSP2_alfa_beta_face_weights` supplies these weights to the
thermodynamic, turbulent, luminosity, and scale-height closures. The shared
`auto_diff_support::get_RSP2_conv_velocity` path now uses the same choice.
Face-to-cell averages and geometric half-cell widths remain centered because
they are not cell-to-face interpolations.

At `k = 1`:

- `RSP2_use_L_eqn_at_surface = .true.` retains
  \(L=4\pi r^2 c a T^4\,\mathtt{RSP2\_Lsurf\_factor}\).
- `RSP2_use_L_eqn_at_surface = .false.` uses the normal atmosphere temperature
  boundary condition in `hydro_eqns::PT_eqns_surf`. The surface radiative
  luminosity is then the model luminosity, `Lr(1) = L(1)`. The existing
  `use_RSP_L_eqn_outer_BC` and `constant_L` controls retain their normal
  precedence when selected.

The explicit RSP2 surface equation takes precedence over `constant_L` at the
surface. Interior `constant_L` rows are unchanged.

The `use_dPrad_dm_form_of_T_gradient_eqn` path uses the reconstructed RSP2
radiative luminosity `Lr_ad`. It therefore remains consistent with the
algebraic flux partition instead of substituting the total luminosity.

## Linear analysis

`star_LNA` now linearizes the same equation selection as the stellar solve:

- Interior RSP2 luminosity slots use the temperature-gradient residual with
  the AD value in `s% gradT_ad`,
  \[
  R_T=\Delta m
  \left(\frac{d\ln P}{dm}\right)_{\rm QHSE}\nabla_T
  -(\ln T_{k-1}-\ln T_k)=0.
  \]
- The surface uses the RSP2 free-streaming row or the atmosphere temperature
  row according to `RSP2_use_L_eqn_at_surface`.
- The `dPrad/dm` temperature row uses `s% Lr_ad` for RSP2.
- RSP2 turbulent-energy and flux audits retain the shared `Y`, `Lr`, `Lc`, and
  `Lt` closures.

## Restart behavior

This reformulation should reduce the thermal discontinuity when changing a TDC
photo to RSP2. The supplied luminosity enters the local gradient inversion, so
the first stellar residual measures temperature-gradient consistency rather
than a large direct flux mismatch.

It does not make the complete TDC-to-RSP2 switch continuous. The following
remain separate:

- RSP2 turbulent velocity initialized from `mlt_vc` is not necessarily on the
  RSP2 turbulent-energy branch.
- Changing `TDC_alpha_M` to `RSP2_alfam` changes the eddy-viscous momentum
  source immediately.
- The RSP2 zero-to-positive turbulent branch still needs a valid branch
  selection rule.
- Remeshing changes the spatial state but does not reconcile the TDC and RSP2
  closures by itself.

## Low-velocity initialization is separate

`RSP2_adjust_vars_before_call_solver` originally treated
`RSP2_w_min_for_damping` as an initialization threshold. With the default
controls, every accepted turbulent velocity below `100 cm/s` could be replaced
by the local steady root before another stellar solve. This conflated damping
regularization with branch initialization.

The discarded replacement repeatedly solved a scalar thermal closure while
changing `w`. It still did not prevent accepted positive-`Y_face`, zero-`w`
cells, and it omitted the coupled turbulent-pressure, turbulent-flux, and
eddy-viscous terms in the stellar row.

The implemented initializer acts only when the accepted state has exactly
`w_start = 0`. For `RSP2_alfat = 0`, the local backward-Euler source and damping
terms factor as

```math
w\left[B\Delta t\,w^2+(1+C\Delta t)w-A\Delta t\right]=0,
```

so the positive initial estimate is

```math
w_{\rm onset}=
\frac{2A\Delta t}
 {1+C\Delta t+
  \sqrt{(1+C\Delta t)^2+4AB\Delta t^2}}.
```

`A`, `B`, and `C` use the same buoyant source, cubic dissipation, and radiative
damping coefficients as `compute_Source`, `compute_D`, and `compute_Dr`.
They are evaluated once from the pre-solver state. Every nonzero accepted
`w_start` is preserved, including a small remap tail. Such a state does not
have the exact zero root and is advanced by the full stellar equation.

The initializer is only an initial guess. The full solve retains
turbulent-pressure work, turbulent-flux divergence, and eddy-viscous work.
`RSP2_w_min_for_damping` has been removed, and both nonlinear RSP2 and
`star_LNA` now use the physical cubic dissipation `B*w**3`. RSP2 includes
turbulent energy independently of `TDC_include_eturb_in_energy_equation`.

## Live-run result after the closure-consistent seed

The 2026-08-29 run in `dev_TDC_delta_scuti` did not retain the positive branch.
Despite the directory name, the active inlist loaded `standard_he_dep.mod` and
evolved a 5.91 solar-mass Cepheid envelope. In profile 6000, faces 73 through
90 had positive `Y_face`, increasing from 0.04 to 35.5, while the RSP2
convective velocity was exactly zero. Face 91 jumped to 3.88 km/s with
`Lc/L = 0.942`. In profile 7000, a smaller zero-velocity, positive-`Y_face`
region had moved to faces 21 through 25.

The surface luminosity excursion near models 6950 through 6960 occurred with
almost no radius change. `log_L` changed by about 0.10 dex while the radius
remained near 73.662 solar radii and the surface velocity remained near
1 km/s. The profile plots showed simultaneous `Eq` and `Uq` spikes. The branch
seed therefore did not solve the accepted zero-velocity state and is not a
complete explanation for the surface thermal excursion.

Before the current correction, the remaining zero branch followed from the
RSP2 turbulent-energy row scaling. `do1_turbulent_energy_eqn` called
`set_energy_eqn_scal`, for which

```math
{\tt scal}=\frac{\Delta t}{e_{\rm start}}.
```

The subsequent division by `dt` makes the dimensionless RSP2 row

```math
R_w=\frac{\Delta w^2+P_t\Delta V+
\Delta t\,\nabla_m L_t-\Delta t(C+E_q)}{e_{\rm start}}.
```

At `w = 0`, the source and eddy-viscous terms vanish because they are
proportional to `w`. If `w_start = 0`, the remaining residual is approximately

```math
R_w(0)\simeq
-\frac{\Delta t\,D_0 w_{\min}^3}{e_{\rm start}},
```

which is small compared with the gas internal energy even when `Y_face` is
strongly positive. The generic stellar solver can therefore accept the
numerical zero branch. Changing the Newton initial guess cannot correct that
convergence test.

At face 90 of profile 6000, `energy = 3.90d12 erg/g`, the timestep was
75.8 seconds, and the pressure scale height was 0.318 solar radii. With the
default `w_min = 100 cm/s`, the scaled damping-floor residual was only
`1.27d-15`, compared with the active `1d-9` maximum residual tolerance.

RSP uses the same steady source and damping estimate in `check_omega`, but its
dedicated turbulent-energy row is assembled without `set_energy_eqn_scal`.
The former RSP2 implementation also used a steady-root initial guess, but its
finite-difference `Y_face` and luminosity equation gave a different nonlinear
path. The closure-consistent backward-Euler seed is useful as a diagnostic, but
it should not be retained as the primary fix. It can also be a poor initial
guess in a shock because it omits `Eq`, which is large in the failing profiles.

The current source gives this row a fixed turbulent-energy scale and retains a
direct onset estimate only for an exactly dormant accepted cell. A new restart
run is required to test the correction. The current run also leaves
`RSP2_num_outermost_cells_forced_nonturbulent = 1` commented out, while RSP
zeros its outer turbulent cell and every established RSP2 test case sets this
control to 1. That boundary difference can affect the photosphere, but it does
not explain the interior zero-velocity strip by itself.

## Discrete radiative equivalence

For the active `dPrad/dm` temperature row,

```math
P_{\rm rad,k-1}-P_{\rm rad,k}
=-\frac{\Delta m_k\kappa_k}{cA_k^2}L_{r,k}.
```

The RSP2 closure supplies `Lr = Lrad_coeff*gradT`. Eliminating `gradT` from
these two relations gives the same finite-difference radiative luminosity as
the former interior RSP2 luminosity row,

```math
L_{r,k}=-\frac{cA_k^2}{\Delta m_k\kappa_k}
\left(P_{\rm rad,k-1}-P_{\rm rad,k}\right).
```

This identity is exact when both paths use the same face opacity and the
`min_kap_for_dPrad_dm_eqn` floor is inactive. The current default
`RSP2_use_mass_interp_face_values = .true.` supplies the same interpolation in
the audited run. The change in nonlinear behavior therefore comes from the
coupled `Y_face`, convective-flux, and turbulent-energy branch, not from losing
the old radiative finite difference.

## Nonlinear Cepheid run audit

The 2026-08-29 run in `dev_TDC_delta_scuti` loads `standard_he_dep.mod`.
This is the 5.91 solar-mass, 75 solar-radius Cepheid envelope, not a delta
Scuti model. The run reached model 20040 and 359.59 days. The first fourteen
measured cycles retained a period near 12 days. Additional radius extrema then
caused the cycle diagnostic to report periods near 6 and 3 days. The maximum
velocity reached `2.446*csound` at model 9314.

The accepted models expose a zero-to-positive RSP2 turbulent branch problem.
In profile 10, at model 9000, face 90 has

```math
\mathcal{Y}=94.54,
\qquad
e_{\rm turb}=0,
\qquad
L_c=0.
```

At face 91 the RSP2 convective velocity is 6.12 km/s and convection carries
94.8 percent of the luminosity. Similar one-zone transitions recur in the
profiles near models 12000, 15000, and 19000. The solver trace from models
18541 through 20040 contains 48 retries, 40 from residual failure, and usually
places the largest correction in `w` near the moving convection front.

The sibling TDC run of the same Cepheid remains on one approximately 12-day
cycle through model 23401. Its dynamical pressure-gradient factor remains near
unity in the saved profiles, while the RSP2 profiles develop negative and
large factors after the first nonlinear disturbance. The RSP2 run still has a
run energy error near `1d-7` and a median of three solver iterations. This is
not a global energy failure. It is a local turbulent-branch and shock problem.

The bounded 100 cm/s Newton seed did not keep an unstable cell on the positive
turbulent branch. Newton returned to a numerically accepted zero-`w` state while
the neighboring cell converged to a velocity of several km/s. The replacement
seed follows the `w` dependence of the active thermal closure and is not capped
at 100 cm/s. A restart run is still required to test whether this removes the
one-zone branch discontinuity.

### Reassessment of the local branch seed

The closure-consistent scalar root is not retained as the proposed fix. It
revisits every cell below `RSP2_w_min_for_damping`, varies `w` while holding the
stellar structure and neighboring amplitudes fixed, and repeatedly evaluates
the algebraic thermal closure. It therefore solves neither the complete RSP2
row nor a fixed-coefficient onset equation. The construction does not establish monotonicity of the combined closure.

The saved run also showed that the extra solve did not select the desired
branch. Profiles through model 60000 continued to contain accepted cells with
positive `Y_face` and `w <= 100 cm/s`. Profile 6000 had 19 such cells and
profile 11000 had 22. The scalar seed changed the initial guess, but the
turbulent-energy residual was still normalized by the gas internal energy and
could accept the low-amplitude state.

For the local pulsation case, `RSP2_alfat = 0`, so a newly unstable cell with
accepted `w_start = 0` has the fixed-coefficient backward-Euler factor

```math
w\left[B\Delta t\,w^2+(1+C\Delta t)w-A\Delta t\right]=0.
```

The positive onset estimate is available directly:

```math
w_{\rm onset}=
\frac{2A\Delta t}
 {1+C\Delta t+
  \sqrt{(1+C\Delta t)^2+4AB\Delta t^2}}.
```

Here `A` is the buoyant source coefficient, `B` is the cubic dissipation
coefficient, and `C` is the radiative damping coefficient. Each coefficient
is evaluated once from the pre-solver state. This estimate is an initial guess
for an exact dormant cell, not a separate closure solve. Every accepted
positive `w_start` must be preserved.

The fixed turbulent-energy scaling experiment was tested and rejected. Before
this branch, the active RSP2 row used `set_energy_eqn_scal`. For an interior
cell its dimensionless residual is

```math
R_{w,k}=
\frac{\Delta e_{{\rm turb},k}+P_{{\rm turb},k}\Delta V_k
      +\Delta t\,\partial_m L_{{\rm turb},k}
      -\Delta t\,C_k-\Delta t\,E_{q,k}}
     {e_{{\rm start},k}}.
```

Here `energy_start` is the gas internal energy. The surface row has the
standard `1d-6` factor, and `dedt_eqn_r_scale` can reduce the scale further.
Old RSP2 used `w/csound = 0` only for cells forced to be nonturbulent. It did
not use `csound**2` to scale an active turbulent-energy row.

Two replacement scales failed in `dev_TDC_RSP2_Cepheid`. Scaling an exactly
dormant row by its first assembled residual normalized that residual to one.
Replacing that fallback with `csound_start**2` did not fix zone 133 because
the accepted `w_start` there was small and positive. Its active-row scale was
therefore `w_start**2`, which again produced an order-unity residual when the
trial solution returned toward zero. The transient scale array has been
removed and the standard MESA energy scaling restored. The analytic onset
value remains only an initial guess for an exactly dormant cell.

### Positive-root selection

The standard scaling did not remove the zero branch. In the accepted model
5000, cell 96 had

```math
Y=8.43, \qquad v_{\rm MLT}=6.74\ {\rm km\ s^{-1}}, \qquad w=0,
```

while cell 97 was on the positive turbulent branch. The abrupt one-cell
transition moved through the envelope and produced the repeated luminosity
and effective-temperature features.

For `RSP2_source_seed = 0`, `RSP2_alfat = 0`, and an accepted
`w_start = 0`, the complete local backward-Euler turbulent-energy residual
factors exactly as

```math
F(w)=wG(w).
```

The terms divided by `w` are formed analytically in
`hydro_rsp2::do1_turbulent_energy_eqn`:

```math
G(w)=w+
\frac{P_{\rm turb}\Delta V}{w}
-\Delta t\left(
\frac{S-D-D_r}{w}+\frac{E_q}{w}
\right).
```

No numerical division by `w` is used. `compute_Source_div_w`,
`compute_D_div_w`, `compute_Dr_div_w`, and `compute_Eq_div_w_cell` retain the
same AD dependence as the corresponding full terms. The linear driving at
the zero solution is

```math
A=\left.\left(\frac{S}{w}+\frac{E_q}{w}\right)\right|_{w=0}.
```

When `A > 0`, the zero root is unstable and the solver uses `G = 0`, which
selects a positive root of the existing equation. When `A <= 0`, the solver
keeps the original `F = 0` row and the radiative solution remains available.
The test is repeated from the current Newton state, so a cell can become
unstable during the solve. The accepted `w_start` is not changed and
`RSP2_source_seed` remains zero.

The branch selection is restricted to `RSP2_alfat = 0`. With nonlocal
turbulent luminosity, the cell-local divergence of `L_t` need not contain the
same factor of `w`, so this deflation is not exact. The change does not alter
the positive-root equations and therefore does not require a separate
`star_LNA` closure change.

The shifted dissipation `B*(w**3 - w_min**3)` injected turbulent energy below
`w_min`. `compute_D` and `rsp2_damping_for_star_LNA` now both use `B*w**3`, and
the obsolete `RSP2_w_min_for_damping` control has been removed.

This is the smallest branch correction supported by the current pulsation
evidence. The archived `asinh`, Fischer--Burmeister, active-set, continuation,
and coupled predictor experiments address exact support during long-term
evolution and are not part of this pulsation fix.

Several active controls also make this run more aggressive than the existing
6 solar-mass RSP2 example. The current soft sound-crossing limit is 4 and the
hard limit is disabled. The existing example uses 0.5 and 4, respectively,
with `max_timestep = 1200`. It also forces the outermost cell nonturbulent and
uses the zero-gas-pressure surface boundary. The current run instead uses the
atmosphere momentum boundary. It also combines Ledoux convection with
`TDC_use_dynamical_gradL`, whose composition treatment remains incomplete.

The smallest discriminating restart sequence is:

1. Restart from `x00008000` and change only the soft and hard sound-crossing
   limits to 0.5 and 4, with `max_timestep = 1200`.
2. If the branch discontinuity remains, repeat with
   `RSP2_num_outermost_cells_forced_nonturbulent = 1` and the established RSP2
   surface pressure boundary.
3. Compare `TDC_use_dynamical_gradL = .false.` before changing the RSP2
   temperature-gradient implementation.
4. Add focused `w`, `Y_face`, `L_c/L`, and pressure-gradient-factor output
   across the moving convection front before changing the branch rule.

## Implementation checklist

- [x] Add a private AD `compute_RSP2_gradT` helper.
- [x] Derive and reuse one `K_r` in the gradient closure and `L_r` diagnostic.
- [x] Derive and reuse one linear `K_c` in the gradient closure, `PII`, and
  `L_c` diagnostic.
- [x] Include the current nonlocal `L_t` in the gradient numerator.
- [x] Honor `TDC_use_dynamical_gradL` through
  `get_TDC_dynamical_gradL`.
- [x] Evaluate the dynamical pressure-gradient factor entirely on the current
  grid so TDC and RSP2 do not mix pre-remesh and post-remesh states.
- [x] Complete a MESA install after the current-grid and radiative-branch
  corrections.
- [x] Avoid subtracting and adding the large Ledoux term on the radiative TDC
  branch.
- [x] Use the resulting neutral gradient in every RSP2 `Y_face` consumer.
- [x] Route interior RSP2 `equL` rows to the MESA temperature equation.
- [x] Select the RSP2 luminosity or atmosphere boundary equation at `k = 1`.
- [x] Use mass interpolation consistently for RSP2 cell-to-face values.
- [x] Remove the alternate RSP definition of `Y_face`.
- [x] Match the RSP2 equation selection in `star_LNA`.
- [x] Test and reject the bounded and closure-consistent scalar seeds.
- [x] Confirm from profiles 6000 and 7000 that the seed does not prevent the
  accepted zero-velocity, positive-`Y_face` branch.
- [x] Remove the closure-consistent scalar root.
- [x] Preserve every nonzero accepted `w_start` and use the analytic
  finite-step onset estimate only for an exact dormant cell.
- [x] Test and reject fixed turbulent-energy scaling for active and dormant
  RSP2 rows.
- [x] Restore `set_energy_eqn_scal` for the active RSP2 turbulent-energy row.
- [x] Select the positive root of the exact local zero-seed equation when its
  zero solution is linearly unstable.
- [x] Keep `w_start` unchanged and retain `RSP2_source_seed = 0`.
- [x] Remove the shifted `w_min` dissipation from RSP2 and make the same change
  in the RSP2 `star_LNA` closure.
- [x] Install MESA after the RSP2 seed correction.
- [x] Reinstall MESA after removing the `w_start` branch test.
- [x] Install MESA after the closure-consistent branch seed on 2026-08-29.
- [x] Install MESA after the final onset and row-scaling correction.
- [x] Reinstall MESA after adding the dormant `csound_start**2` scale on
  2026-08-30.
- [x] Reinstall MESA after restoring the standard energy-row scaling on
  2026-08-30.
- [ ] Run an RSP2 restart after adding the positive-root equation row.
- [ ] Confirm that positive-`Y_face`, zero-`w` cells and the moving one-cell
  luminosity feature are absent.
- [ ] Check AD partials for `lnT`, `lnd`, `lnR`, `L`, `w`, and neighboring
  variables.
- [ ] Compare direct and reconstructed `L_r + L_c + L_t` at convergence.
- [ ] Test a static RSP2 model, a TDC-to-RSP2 photo switch, a pulsating
  restart, `TDC_use_dynamical_gradL` on and off, Schwarzschild and Ledoux
  cases, and remesh on and off.
- [ ] Form and smooth a separate `B_QHSE` before claiming complete Ledoux
  support for the dynamical neutral gradient.

## Source references

- `star/private/hydro_rsp2.f90`
  - `do1_rsp2_L_eqn`
  - `do1_turbulent_energy_eqn`
  - `compute_Y_face`
  - `compute_PII_face`
  - `compute_D`
  - `compute_Lr`
  - `compute_Lc`
  - `compute_Lt`
  - `RSP2_adjust_vars_before_call_solver`
- `star/private/star_LNA_turbulence_closures.f90`
  - `rsp2_damping_for_star_LNA`
- `star/private/hydro_eqns.f90`
  - interior RSP2 `equL` routing
- `star/private/hydro_temperature.f90`
  - `do1_dlnT_dm_eqn`
  - `do1_alt_dlnT_dm_eqn`
  - `eval_dlnPdm_qhse`
- `star/private/auto_diff_support.f90`
  - `get_RSP2_conv_velocity`
- `star/private/star_LNA_support.f90`
  - `assemble_luminosity_rows`
  - `temperature_gradient_resid_for_star_LNA`
  - `dPrad_dm_resid_for_star_LNA`
- `star/private/turb_support.f90`
  - `get_TDC_dynamical_gradL`
- `turb/public/turb.f90`
  - `set_TDC`
- `star/defaults/controls_dev.defaults`
  - RSP2 closure and initialization controls

The earlier Flaskamp development record is
`notes/EbF/rsp2/old_notes/rsp2_flaskamp_gradient_reformulation.md` in the
clean MESA checkout.
