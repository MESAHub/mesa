# RSP2 local three-equation implementation plan and work record

Started 2026-09-19 on `EbF/star_lna`, base `4b5e4e999`.

**Current status: the cell-w implementation is checkpointed in `d28e05018`.**
The subsequent face-w implementation for RSP2 and RSP3 is tracked in
[rsp2_face_w_implementation.md](rsp2_face_w_implementation.md). It remains
uncommitted, is compiled and installed, and has its first stellar test recorded there. That record supersedes the
cell/face placement and remapping descriptions below.

The cell-w run history below includes an unresolved variance closure failure.
The separate [covariance correction](rsp3_covariance_closure.md) remains
unimplemented. Moving w to faces changes placement, not that closure.
Prior installation and run records apply to the cell-w checkpoint. The shared
gradient helpers were committed as `1700f15ac` before the optional mode.
No user case was relinked and no push has occurred.

The remaining record documents the earlier implementation and its checks.

## 1. Scope and invariants

Add `RSP2_use_3equation_model = .false.` to `&controls`. This is a persistent
choice of convection equations, not a one-time star_job action. Existing
star_job controls still enable RSP2 itself. Turning the new control on while
RSP2 is off does not silently turn on hydrodynamics.

The two new moments and their evolution equations live on FACES:

- `Pi = <v_r' s'>`, signed.
- `Phi = <s'^2>`, the full nonnegative variance.

Cell `w` with `e_t=w^2`, face Y, and face L retain their existing locations.
`k=1` is the surface; k increases inward. `PII_face` remains the old factored
one-equation quantity, with different dimensions from the new entropy flux.
No extra factor of w multiplies Pi in luminosity. The cell buoyancy source
uses the dimensionless reconstruction factor w_cell/sqrt(e_t_face) below.

| Quantity | Location | cgs units | Role |
| --- | --- | --- | --- |
| `w`, with `e_t=w^2` | Cell | cm/s | Existing turbulent energy state; zero state branches are specified below. |
| `Pi = <v_r' s'>` | Face | erg cm/(g K s) | New signed solver variable; supplies Lc directly. |
| `Phi = <s'^2>` | Face | erg^2/(g^2 K^2) | New nonnegative solver variable; full variance, not half-variance. |
| `Y_face` | Face | Dimensionless | Existing solver variable in `gradT=gradL+Y_face`. |
| `L` | Face | erg/s | Existing total luminosity variable and balance row. |

"Three equations" means three turbulence evolution equations. The existing
algebraic flux-balance and temperature-gradient rows remain additional rows.

Preserve when the control is false:

1. Existing RSP2 variable/equation counts and order.
2. Existing PII(Y), source, Dr, Lt, Eq/Uq and physical energy/gradient equations.
3. TDC equations and TDC face velocity placement.
4. Existing old model/photo readability, including one-equation RSP2 photos.
5. User inlists, cases and outputs; do not toggle their new control for them.

Do not add a physical enthalpy-flux cap, alter luminosity time centering,
change Eq/Uq, or fix the unrelated k+2 conservative-work Jacobian in this diff.
The subsequent request to address zero-w handling in BOTH modes is recorded
in section 2.4. It permits a shared numerical residual/domain treatment;
it does not authorize changing those physical terms or time weights.

## 2. Equations and numerical choices

The literature derivation, normalization and references are in
[the local-moment proposal](rsp2_three_equation_local_proposal.md). That remains the theoretical
record; this file tracks implementation against the actual MESA call paths.

In uniform composition the continuous thermal-moment equations are

```math
D_t Pi=(2/3)e_t(-d_r s)
       +[-(d_r P)/rho](chi_T/chi_rho)/c_p <s'^2>
       -(alfa_pi sqrt(e_t)/Lambda+1/tau_rad)Pi
       -Pi d_r u,
```

```math
D_t <s'^2>=2(-d_r s)Pi
          -(alfa_phi sqrt(e_t)/Lambda+2/tau_rad)<s'^2>.
```

The final mean-strain term in the flux equation is present in the full RANS
moment equation and absent from the quiet-background Flaskamp reduction.
The code retains this term with the u face difference and interpolated v
cell slopes specified in the implementation record. This does not introduce a new viscosity
coefficient. No nonlocal Pi or variance transport is included. Lt transport
in the kinetic-energy equation stays active for nonzero RSP2_alfat.

```math
1/tau_rad = 4 sigma gamma_r^2 T^3/(c_p kappa rho^2 Lambda^2),
L_c = 4 pi r^2 (rho T)_f Pi_f,
C_k = Source_k - D_k,
D_k = RSP2_alfad (8/3)sqrt(2/3) w_k^3/Lambda_k.
```

The direct one-equation `Dr=e_t/tau_rad` is disabled only in the new mode.
Radiative damping appears in the two entropy moments instead. No additional
radiative sink may be silently added to the total energy equation.

New coefficients:

- `RSP2_alfa_pi = 1d0` multiplies `x_ALFAPI = 6*sqrt(2/3)`.
- `RSP2_alfa_phi = 1d0` multiplies `x_ALFAPHI = 4*sqrt(2/3)`.

The alpha symbols in the equations are the products of the control and source
constant. The constants are the local MLT calibration listed in Braun et al.
(2026), section 2. Her calibrated solar models use different coefficients;
unit controls do not reproduce that calibration. At the user's request,
only these two new defaults are changed. Existing RSP2 defaults are retained.

Use the existing gamma_r and the public TDC cell/face scale-height and
mixing-length helpers. Do not introduce another Hp solver variable or
another set of scale-height controls.

The controls documentation also gives the normalization relative to Braun
et al. (2026), equations 9 and 10 and Table 1:

| Control | Meaning here | Teresa's value in this normalization | Current default |
| --- | --- | --- | --- |
| `RSP2_alfar` | `gamma_r = 2*sqrt(3)*RSP2_alfar`; cooling rate is proportional to its square | `1d0`, giving gamma_r about 3.46 | `0d0` |
| `RSP2_alfat` | coefficient of `Lambda*sqrt(e_t)` in kinetic energy diffusion | `0.25d0`, equal to alpha_omega with no extra velocity factor | `0d0` |

In the one equation model alfar gives cell `Dr=e_t/tau_rad`. In the three
equation model direct Dr is absent and the face moments receive
`-Pi/tau_rad` and `-2*Phi/tau_rad`. Alfat transports kinetic energy in both
modes, without adding moment transport. These coefficient translations do
not assert identical finite mesh interpolation or reproduce her full model.
The user requested explanatory comments, not changes to these two defaults.

### 2.1 Entropy-gradient audit: preserve the discrete neutral state

**Previous logarithmic choice withdrawn. The replacement below is implemented.** The
user identified the risk of noisy driving during evolution. Checking the
neutral state confirms a stronger issue: exact inversion of the temperature
row followed by a logarithmic entropy difference introduces a finite-zone
bias even with smooth data, exact hydrostatic balance and high-precision
arithmetic. This is not cured by `log1p` or tighter Newton tolerances.

For the standard temperature row, take uniform composition, constant
grad_ad, equal cell masses, no reconstruction/time centering, and pressure
values satisfying the discrete hydrostatic balance. Then

```math
(P_k-P_{k-1})/Ppoint = 2 tanh[ln(P_k/P_{k-1})/2],
(T_k-T_{k-1})/Tpoint = grad_ad (P_k-P_{k-1})/Ppoint,  Y_face=0.
```

The logarithmic inversion below infers the apparent superadiabaticity

```math
ln(T_k/T_{k-1})/ln(P_k/P_{k-1}) - grad_ad
= grad_ad (grad_ad^2-1) [ln(P_k/P_{k-1})]^2/12
  + O([ln(P_k/P_{k-1})]^4).
```

For grad_ad=0.4 and ln(P_k/P_{k-1})=0.01 this is -2.80e-6
at Y=0. With Y=1e-8 the inferred value is still negative, -2.79e-6.
These are manufactured examples, not values measured in the user's model.
Unequal cell weights can introduce a leading error proportional to the
pressure log jump. Smooth bias can become mesh-dependent variations as
zone sizes and weights change; the check does not demonstrate random noise
or an actual MESA instability. The risk depends on the size of the physical
driving, not simply on whether the model evolves or pulsates.

The replacement discretizes the thermodynamic differential in the selected
temperature row's coordinate and cancels its thermal neutral reference
analytically before adding Y. It changes the finite-zone entropy quadrature;
it is not an algebraically identical rewrite of the rejected logarithmic
gradient. No smoothing, clipping or evolution-specific closure is added.

#### Selected mapping for the optional three-equation mode

Use a consistent face EOS state (section 2.2) and the selected temperature
row to compute `temperature_gradient_per_gradT`, with units 1/cm. This is
a coefficient, not a solver variable. The row implies
`-dlnT/dr = temperature_gradient_per_gradT*(gradL+Y_face)` at convergence,
where the left side denotes the row's discrete thermodynamic differential:
`(T_k-T_(k-1))/T_face` for standard, the log ratio for actual-log, and
`(Prad_k-Prad_(k-1))/(4*Prad_face)` for dPrad, each multiplied by
`4*pi*r^2*rho_face/dm_temperature`. These are different finite-zone
quadratures of the same continuous derivative.
The physical thermal pressure gradient uses CURRENT Peos, not turbulent
pressure: `pressure_gradient = -dln(Peos)/dr` in the chosen face quadrature.

For standard and radiation-pressure rows use

```math
pressure_gradient = (4*pi*r^2*rho_face/dm_temperature)
                    *(Peos_k-Peos_(k-1))/Peos_face.
```

For the actual-log row use the pressure log ratio with the same geometric
factor instead. The temperature coefficient is:

| Temperature row | temperature_gradient_per_gradT |
| --- | --- |
| Standard | `4*pi*r^2*rho_face*(-dlnPdm_qhse)*(Tpoint/T_face)` |
| Actual log | `(4*pi*r^2*rho_face/dm_temperature)*ln(Peos_k/Peos_(k-1))` |
| dPrad | `rho_face*kap_row*Lrad_per_gradT/(4*clight*area*lambda*Prad_face)` |

`Tpoint` is the temperature row's mass-weighted T. `Peos_face`, T_face and
Cp_face are the common moment/EOS face state. Retain the ratios when the
reconstruction or RSP2 face weights differ from the temperature row's
weights. The dPrad coefficient uses the SAME opacity floor, radiation flux
factor and inner-boundary spacing as that row; Prad_face=crad*T_face^4/3.
The common RSP2 radiative coefficient in the new mode uses that face state.

For standard/dPrad, `pressure_gradient_reference` is the actual thermal
pressure gradient when dynamical gradL is on. When off, use the expected
HSE pressure jump from the existing gravity/area operator, with the
momentum face mass correction when enabled, divided by the temperature
stencil mass and Peos_face. Thus time weights, rotation, mass correction
and the Tpoint/Ppoint ratios are explicit. In the actual-log row the
reference is always its actual pressure gradient: the coordinate already
uses the resolved pressure difference, so no HSE conversion is required.

Set the new mode's neutral coordinate and thermal entropy driving together:

```math
gradL = (grad_ad_face + gradL_composition_term)
        *pressure_gradient_reference/temperature_gradient_per_gradT,

thermal_entropy_gradient = Cp_face *[
   temperature_gradient_per_gradT*Y_face
   + pressure_gradient_reference*gradL_composition_term
   + grad_ad_face*(pressure_gradient_reference-pressure_gradient)].
```

When the reference is the actual pressure gradient, cancel the last term
ANALYTICALLY. Do not form gradT-gradL or subtract the large adiabatic
terms. For the actual-log row, set gradL=grad_ad_face+composition directly;
this also avoids 0/0 when the pressure jump vanishes. For standard/dPrad
an invalid temperature coefficient is an error, not a denominator floor.
When dynamical gradL is off, retain the physical HSE-versus-actual pressure
term. Its numerical precision is limited by the resolved pressure force;
it must not be silently dropped to manufacture exact neutrality.

This is an explicit change of Y's neutral coordinate in THREE-EQUATION
mode. Ordinary RSP2 and TDC retain their current gradL path. At fixed
physical gradT, changing the new mode's reference and shifting Y by the
opposite amount leaves the thermal driving unchanged. On mode conversion,
initialize `Y_new=gradT_old-gradL_new` to preserve the old gradient before
the coupled solve; preserve histories under the same coordinate convention.

With varying composition this is the THERMAL part of the entropy gradient.
Use the existing `gradL_composition_term` from `hydro_vars:set_grads`,
including its existing Brunt-B policy; add no new smoothing. Its conversion
is not a model for composition fluctuations or a full Ledoux moment closure.
Existing real-valued composition coefficients remain frozen with respect to
the thermal AD variables, as in the present infrastructure.

`constant_L` replaces the temperature row. In that path evaluate the measured
thermal differential from cell T and Peos with the common face denominators;
Y is not a surrogate for the missing temperature equation. Its tiny-driving
accuracy is limited by the cell/EOS values. Surface and forced moment rows
retain their explicit zero-moment boundary conditions.

Implementation placement: `star/private/hydro_gradient_support.f90` shares
QHSE/gravity and radiation-row factors without the dependency cycle
`hydro_temperature -> hydro_momentum -> hydro_rsp2 -> hydro_temperature`.
Existing QHSE/gravity behavior is the default; the optional instantaneous
path is now used by LNA. Source call sites are:

- `hydro_rsp2:compute_RSP2_gradT`: new-mode gradL and the existing Y coordinate.
- `hydro_rsp2:rsp2_moment_rhs`: the same thermal driving in both rows.
- `hydro_rsp2:compute_Lrad_coeff`: common face EOS/opacity in the new mode.
- `hydro_temperature:do1_alt_dlnT_dm_eqn`: shared opacity floor, radiation
  factor and boundary spacing, with the existing residual unchanged.

Ordinary RSP2/TDC keep their gradL and radiative coefficient paths. The new
mode uses the shared mapper after conversion of Y, including in its LNA rows.

#### Previous logarithmic candidate retained for the derivation record

The rejected choice used the thermal entropy change between the two cells,
with uniform-composition thermodynamics:

```math
(-d_r s)_f = c_{p,f} (4 pi r_f^2 rho_f/dm_bar_f)
 [ln(T_k/T_{k-1}) - grad_ad,f ln(P_k/P_{k-1})].
```

Here k is the inner cell, k-1 the outer cell. Positive driving means entropy
decreases outward. This is a two-point discretization of the thermodynamic
differential, not an exact integration through an arbitrarily varying EOS.

Its numerical evaluation replaces the temperature jump by the expression
implied by the selected temperature residual. This gives the same coupled
roots and, at a root, an equivalent linearized system after a row operation.
It avoids subtracting two large cell entropies, but still combines different
finite-zone quadratures. Row equivalence proves neither neutral-state
preservation nor accuracy relative to a tiny physical driving term. The
replacement is a Newton residual choice, not a new physical source. Away
from a root the residual/Jacobian can differ.

Use `log1p((P_k-P_{k-1})/P_{k-1})` for the pressure log ratio, sharing the
same current EOS pressure jump with buoyancy. This avoids subtracting large
absolute logarithms; it cannot recover precision already lost in EOS values.
`auto_diff_real_star_order1` already supplies `log1p` and `expm1` with AD.

**Actual-log-gradient temperature row.** For
`use_gradT_actual_vs_gradT_MLT_for_T_gradient_eqn`, the thermal difference is

```math
ln(P_k/P_{k-1}) [Y_face + (gradL-grad_ad,f)].
```

Evaluate it in this order. Do not first form `gradT=gradL+Y_face` and then
subtract grad_ad: that can erase Y in floating point. The existing gradL
may include its dynamical/composition factors; do not apply those again.

**Standard temperature row.** The code imposes

```math
(T_k-T_{k-1})/Tpoint = -dm_bar*dlnPdm_qhse*(gradL+Y_face),
Tpoint = a*T_k+(1-a)*T_{k-1},
a = dm(k-1)/(dm(k-1)+dm(k)).
```

Here a is the actual temperature-row weight, regardless of the RSP2 face
weight option. For clarity define `temperature_jump` as this dimensionless
ratio; it is a computed quantity, not an additional solver variable:

```math
ln(T_k/T_{k-1}) =
 log1p((1-a)*temperature_jump)-log1p(-a*temperature_jump).
```

To evaluate the small thermal departure accurately, compute the jump for
an adiabatic temperature ratio, then its difference from the predicted jump:

```math
temperature_jump_ad =
 expm1(grad_ad,f*ln(P_k/P_{k-1})) /
 [1+a*expm1(grad_ad,f*ln(P_k/P_{k-1}))],

delta_temperature_jump = -dm_bar*dlnPdm_qhse*Y_face
                       +(-dm_bar*dlnPdm_qhse*gradL-temperature_jump_ad).
```

The bracket in the entropy-gradient formula is then

```math
log1p((1-a)*delta_temperature_jump/[1+(1-a)*temperature_jump_ad])
-log1p(-a*delta_temperature_jump/[1-a*temperature_jump_ad]).
```

Reuse `hydro_temperature:eval_dlnPdm_qhse`, including its actual pressure,
rotation, gravity and time-centering choices. Do not substitute g/Hp or
copy only its hydrostatic special case. Reject invalid logarithm domains as
invalid trial states, rather than clipping the thermal gradient.

**Radiation-pressure temperature row.** The code imposes

```math
(T_k^4-T_{k-1}^4)/T_{k-1}^4 =
 dm_bar*kap_face*Lr/[clight*(4*pi*r_f^2)^2*flxLambda*Prad_(k-1)].
```

Use the same opacity floor, optional radiation flux factor and inner-boundary
spacing as `do1_alt_dlnT_dm_eqn`. When radiation flux limiting is disabled,
flxLambda is one. This is the existing RADIATIVE limiter, not an enthalpy cap.
`Lr=compute_Lrad_coeff*(gradL+Y_face)` in the current RSP2 code.
Compute the predicted relative radiation-pressure jump minus
`expm1(4*grad_ad,f*ln(P_k/P_{k-1}))`, splitting off the term proportional to Y
before adding gradL. The thermal bracket is exactly

```math
(1/4)*log1p(delta_relative_Pradiation /
                  exp(4*grad_ad,f*ln(P_k/P_{k-1}))).
```

This is an exact inversion of that discrete temperature row, not a continuum
approximation `Delta T4/(4*T_face^4)`. Different gradient forms need not give
identical coarse-grid solutions; each driving term must match its own row.
If a constant-L row replaces the temperature constraint, there is no such
identity to substitute: use the measured thermal log difference in that path.

Possible implementation placement after resolving the discretization: a new entropy-gradient helper in
`hydro_temperature` can reuse its existing QHSE routine and shared radiation
coefficients. Call it from `hydro_eqns` and `star_LNA_support`, then pass its
AD result to the shared `hydro_rsp2` moment rhs. Do NOT make hydro_rsp2 use
hydro_temperature: hydro_temperature already depends on hydro_momentum,
which depends on hydro_rsp2, so that would create a module cycle. Factor out
only radiation-row coefficient work that is genuinely shared by the two
callers; preserve the existing temperature residuals.

A thermal entropy model still omits composition fluctuations. Keeping the
existing Ledoux coordinate does not add their missing moment equations.

### 2.2 Face reconstruction and buoyant source decision

**Implemented with common face thermodynamics.** The thermal mapper, moment
damping, Lrad, Lc and buoyancy use one face state in the new mode:

- With `use_face_reconstruction`, use the existing
  `get_reconstructed_face_eos_kap_ad` result, including Cp, expansion
  coefficients, grad_ad and opacity, with its AD derivatives.
- Otherwise use the existing RSP2 alfa/beta weights on the two cell values.
  Do not use an unrelated helper's unconditional mass weights when
  `RSP2_use_mass_interp_face_values` is false.
- New-mode Lc uses the product `rho_face*T_face`, not the average of two
  cell rho*T products. The old mode retains its existing product average.
  Convert a snapshot with its ACTUAL saved old Lc divided by the new
  `4*pi*r^2*rho_face*T_face`; `average(w)*PII` is not generally the same seed
  after changing this coefficient.
- Reconstruct face kinetic energy as `alfa*w_k^2+beta*w_(k-1)^2`.
  Do not reconstruct it as `square(alfa*w_k+beta*w_(k-1))`.
- Both new moments are already face variables; do not interpolate them
  again inside their own local equations.

Use the pressure force actually resolved between the cells in buoyancy:

```math
buoyancy_f = [4*pi*r_f^2/dm_bar_f]
            [chiT_f/(chiRho_f*Cp_f)] (P_k-P_{k-1}),
Source_k = (w_k/2)[buoyancy_k*Pi_k/sqrt((e_t)_face,k)
                 +buoyancy_(k+1)*Pi_(k+1)/sqrt((e_t)_face,k+1)].
```

`buoyancy_f` is a named physical coefficient, not a new unknown. This avoids
approximating Delta P by P_face*Delta ln P. In hydrostatic conditions it
reduces to g*(chiT/chiRho)/Cp. It does not reinterpret a limited/mixed Hp as
an actual pressure gradient. Keep the public TDC Hp/Lambda functions for
scale heights, turnover and radiative cooling. Reconstruct the face entropy
correlation `Pi_face/sqrt(e_t_face)` to the cell, then multiply by its own w.
This explicitly REPLACES the provisional direct average of buoyancy*Pi:
that direct average can remove energy from an empty cell next to a turbulent
cell when Pi is negative. The revised source vanishes with the cell's w.
It changes the source discretization, not the independently solved face
luminosity or the continuous homogeneous limit. It is a chosen discretization,
not a formula prescribed by the literature's continuous equations.

For ordinary evaluation form the bounded ratio `w_cell/sqrt(e_t_face)` first,
then multiply by Pi, avoiding a large intermediate Pi/sqrt(e_t_face).
The already named source_div_w is needed only for the canceled local row
below. Use the selected exactly dormant branch when both face energies and
its moments vanish; do not manufacture a finite Pi/w from a denominator floor.
For constant buoyancy the reconstructed cell covariance obeys its covariance
bound, using the mean bounding-face variance, IF the face moments obey theirs:
Cauchy--Schwarz gives `(mean(Pi_face/w_face))^2 <= mean((Pi_face/w_face)^2)`.
This interpolation property does not make the evolution closure realizable.
Adjacent-face averaging uses the existing half-cell quadrature; boundary
masks apply before averaging.

The same C=Source-D enters the turbulent-energy equation; both existing gas/
total energy forms already include d(e_t)/dt. No second gas -C is added.
Only two adjacent thermodynamic cells enter each buoyancy coefficient, so
shifting the inner face to the cell equation fits the three-zone AD stencil.
The finalized face-state helper must be used in new-mode grad_ad/gradL too.
Do not change old-mode Eq/Uq or silently widen their stencil.

#### Comparison with current one-equation RSP2

Source audit on 2026-09-19, checked against HEAD as well as the inactive
provisional new-mode branch: current `hydro_rsp2:compute_Source` already
uses the cell's own w. For an interior cell it evaluates

```math
Source_k = (w_k + RSP2_source_seed)
 [P_k*chiT_k/(rho_k*chiRho_k*Cp_k)]
 (1/2)[PII_k/Hp_k + PII_(k+1)/Hp_(k+1)],

PII_f = x_ALFAS*(Lambda_f/Hp_f)*Cp_f*Y_face,f.
```

These are `compute_Source`, `compute_Source_div_w` and
`compute_PII_from_Y` in `star/private/hydro_rsp2.f90`. With the default
`RSP2_source_seed=0`, this buoyant source is exactly zero at cell w=0,
including when a neighboring face has negative PII. It does not average
the complete face heat flux into the cell source. `compute_Lc_terms`
separately uses the interpolated face w, so a shared-face luminosity can
remain nonzero when only one adjacent cell has w=0.

Thus the zero-cell protection motivating the proposed reconstruction is
already present in one-equation RSP2. The proposed independent-moment source
preserves that structural property; it is not a newly discovered missing
factor in the existing source. Its pressure-force coefficient and choice
sqrt(average(w^2)) differ from the old coefficient and average(w), so do
not claim an exact discretization identity between the two modes. A nonzero
source seed deliberately removes the zero-source property and needs its
own physical/branch interpretation.

### 2.3 Normalization decision

**Implemented; runtime tolerance calibration is not yet tested.** Keep
backward-Euler moment source updates at the current COUPL time level:

```math
R_Pi = [Pi-Pi_start-dt*rhs_Pi]/Pi_scale,
R_s2 = [<s'^2>-<s'^2>_start-dt*rhs_s2]/variance_scale.
```

Do not change L/Lt/work time-centering controls. LNA uses continuous rhs and
moment inertia, not these finite-step residuals or their start-state scales.

Tie the flux scale to the existing TDC/RSP2 luminosity scale:

```math
L_scale,f = max(1 erg/s, abs(L_start,f), 1d-3*max_j abs(L_start,j)),
Pi_reference,f = L_scale,f/(4*pi*r_start,f^2*rho_start,f*T_start,f),
Pi_scale,f = max(abs(Pi_start,f), Pi_reference,f).
```

Multiplying Pi_reference by the face conversion coefficient gives L_scale
exactly. This avoids an arbitrary entropy-flux floor in cgs units. For the
variance scale define dimensional reference values, not physical caps:

```math
velocity_reference,f = [L_scale,f/(4*pi*r_start,f^2*rho_start,f)]^(1/3),
variance_reference,f = (3/2)*(Pi_reference,f/velocity_reference,f)^2,
variance_scale,f = max(<s'^2>_start,f, variance_reference,f).
```

The reference speed is the speed with kinetic-flux scale rho*v^3 equal to
L_scale/area. The factor 3/2 converts to the isotropic radial variance. This
is a dimensional normalization, NOT an assertion about the actual velocity
or covariance and NOT an enforced covariance bound. It remains finite at
zero turbulent energy without dividing by the actual w. The final root does
not depend on it; which small residuals satisfy a finite tolerance does.

Compute and store both scales once from the accepted/start face state, after
mode conversion/remesh and before Newton. Use them for the new column scales
as well; normalize corrections by max(start scale, current moment magnitude).
Do not recompute the row scales as the trial moments change. The two scale
work arrays must join alloc's do1 handling and be rebuilt on retry/start;
they are derived state, not extra photo/model columns.

### 2.4 Zero-w decision: retain w and distinguish the physical branches

**Implemented in both RSP2 modes; standalone checks pass. Runtime validation
remains for the user.** Retain cell w so existing Eq/Uq and Lt remain smooth
polynomials/products of the physical turbulent velocity. Do not globally
replace w by energy. That substitution gives a regular local inertia at zero
but produces singular square-root derivatives in the existing terms:

```math
Eq = (Eq_div_w)*sqrt(e_t),
Lt_f = -coefficient_f*[alfa*sqrt(e_k)+beta*sqrt(e_(k-1))]
                        *(e_(k-1)-e_k).
```

The Lt derivative diverges at an interface with only one zero-energy cell;
the Eq derivative diverges with finite mean strain. The u-flag face stresses
have the same issue. Flaskamp section 4.3.2 explicitly discusses why a
sqrt(energy) variable is useful. Substituting zero for infinite derivatives
would not be an exact Jacobian.

These coordinate issues also apply to current one-equation RSP2. The
inertial derivative is `d(w^2)/dw=2*w`, but this alone does NOT imply a
singular full row: the derivative of its source, which is linear in w,
can remain finite at zero. In the local, zero-start, zero-seed case with
no Lt and Eq proportional to cell w, both w=0 and a nonzero solution may
be roots of the energy residual. At zero,

```math
dR_energy/dw = -dt*(Source_div_w + Eq_div_w).
```

If this driving also vanishes, the local row can lose its derivative.
If it is positive, Newton can instead remain on the exact zero root when
the intended solution is the positive branch. These are distinct concerns.

The existing `do1_turbulent_energy_eqn` already cancels the common w factor
when `w_start=0`, `RSP2_source_seed=0`, `RSP2_alfat=0`, the velocity is v
or `RSP2_alfam=0`, and `Source_div_w+Eq_div_w>0`. Its
`RSP2_adjust_vars_before_call_solver` also constructs initial guesses for
positive local buoyancy and incoming Lt. The Lt predictor retains the full
energy residual and does not change accepted/start state. These are limited
existing treatments, not a general zero-state or u-flag heating solution.
No change to them is made by this audit; the new three-equation branch
strategy below is a proposal and has not been installed into either mode.

The source reconstruction in section 2.2 restores a true cell-w factor.
This is why the original statement "the factored row must always be disabled"
is superseded: it applied to the earlier, now-rejected direct face-source
average. With the new reconstruction there are three distinct cases.

**1. Accepted w=0, no Lt, and Eq proportional to that cell's w.** The unscaled
energy residual factors exactly as

```math
R_energy = w*[w + (Ptrb*dV)/w - dt*source_div_w
                  + dt*D_div_w - dt*Eq_div_w].
```

There is no old Dr in this mode. All start-pressure terms proportional to
w_start^2 also vanish here. Evaluate the bracket analytically with the
existing divided routines, never by dividing a computed zero residual by w.
The positive branch sets this bracket to zero. The dormant branch sets w=0.
A concise branch-selection residual is

```math
R_w = min(w, energy_residual_after_canceling_one_factor_of_w) = 0.
```

Both arguments have velocity units; apply the existing velocity/energy-row
normalization afterwards. This is a complementarity choice: w and the
canceled expression are nonnegative, and at least one vanishes. For positive
driving it selects the nonzero branch; for negative driving it permits exact
zero. It removes the double-zero derivative when all driving vanishes.
At a tie use the dormant derivative. This is a piecewise differentiable
active-set residual, not a claim that one smooth Jacobian exists at the tie.
The original energy residual must ALSO satisfy its tolerance at convergence.
It remains exactly zero at either selected root.

The cancellation applies to v-flag cell Eq, and to u-flag when there is no
eddy-viscous term. Do not apply it to u-flag face heating with finite energy
imported from a neighboring cell, or to general nonzero Lt.

**2. Incoming Lt or u-flag viscous production can create cell energy.** Keep
the full existing energy equation, including its time weights. Extend the
existing positive-w initial-guess construction to include the actual imported
Eq and the new moment source as well as Lt. If zero Pi is driven by finite
face energy or variance, predict the face moments first to obtain a usable
initial source; preserve their accepted start values. Freeze coefficients ONLY while
constructing that initial guess. For positive available energy the scalar
quadratic predictor has the form

```math
work_coefficient*w_guess^2-linear_source_coefficient*w_guess
  -available_energy = 0,
```

using the rationalized positive root. Then recompute all current coefficients
and AD derivatives and solve the full residual, including cubic dissipation.
The predictor does not alter w_start, accepted moments, conserved energy,
timestep weights or the final equation. Do not assume that a nonzero alpha_t
turns all source terms into the local factored case. Negative available
energy is not repaired by clipping the equation; use the coupled solve or
reject an inadmissible trial/timestep.

**3. Exactly dormant local region.** With zero accepted moments, zero face
kinetic energy and no current imported energy or shear production, use the
zero-state constraints. Zero face kinetic energy means BOTH adjacent cells
are zero. The face moments' homogeneous evolution from zero has that exact
zero solution, so Pi=variance=0 selects this branch without adding a floor.
Reevaluate the branch as neighboring CURRENT turbulence or physical forcing
changes; do not freeze a whole radiative region for the timestep merely
because its accepted state was zero. Faces adjacent to nonzero energy retain
their moment evolution equations, even when one cell has w=0.

If variance is nonzero at zero kinetic energy, the buoyancy term can generate
Pi: this is NOT a dormant state. A coupled moment/energy initial guess must
include that forcing. Likewise nonzero Pi at zero face energy cannot be
hidden by a zero denominator convention. These are explicit boundary cases
for the prototype and model/remesh validation.

LNA must use the same chosen dormant branch as the nonlinear formulation:
fully dormant face moments get algebraic zero perturbations; mixed faces
retain their dynamical rows. This analysis does not describe spontaneous
activation of an exactly unperturbed dormant region. Preserve the existing
nonzero snapshot turbulence on activation. The old RSP2_source_seed must not
silently change meaning in the new model.

The two-cell prototype below tests these branches with both cell/face
viscous-heating placements, nonzero Lt, positive and negative Pi, an empty
cell next to a turbulent one, and variance remaining at zero kinetic energy.
It solves the coupled two-cell energy and one-face moment system and checks
the original energy/moment residuals after branch selection. It does NOT
solve mean momentum, EOS, radiation, geometry or composition, or exercise
MESA's Newton driver. Those coupled derivatives and code paths remain
implementation/runtime checks; do not infer full MESA convergence from it.

The published reduced isotropic closure still does not guarantee
`Pi^2 <= (2/3)*e_t*<s'^2>`. The proposal records a continuous counterexample.
Nonzero alfat changes kinetic energy without transporting the new moments.
No solver-variable choice, source interpolation or variance clipping proves
realizability or cures a physical branch with no admissible solution.

#### Shared positive-w residual and domain treatment, 2026-09-19

Implemented for BOTH modes: solve the divided physical energy residual at
positive subsonic w. The small-w correction below supersedes the original
implementation that divided the assembled residual. Cancel known local w
factors term by term before differentiation. Lt_theta and Ptrb_theta retain
their EXISTING time weights:

```math
R_energy,k = w_k^2-w_start,k^2 + Ptrb_theta,k*(1/rho_k-1/rho_start,k)
 +dt*[(Lt_theta,k-Lt_theta,k+1)/dm_k-Source_k+D_k+Dr_k-Eq_k],

R_w,k = R_energy,k/w_k,   w_k>0.
```

Dr is present in one-equation RSP2 and absent in the new mode. Division at
positive w is valid with nonzero alfat and u-flag face heating; it does NOT
assert that those terms share a factor w at zero. Mathematically the
derivative is

```math
dR_w/dx_j = (dR_energy/dx_j)/w-(R_energy/w^2)*dw/dx_j.
```

This is more than a fixed row scale away from convergence. At a positive
root it is nonsingular row scaling, so the positive roots are unchanged.
The storage term becomes `w-w_start^2/w`, with derivative
`1+(w_start/w)^2`. Positive incoming energy independent of local w gives
`-incoming_energy/w`, also with positive derivative. This removes the
vanishing storage derivative and helps the known incoming-Lt Newton trap.
It is not a proof of monotonicity or convergence of the full coupled system.

The implemented normalized row equals `R_energy*max(1,csound/w)*scal/dt`.
For `0<w<csound`, it is evaluated as
`[R_energy/w]_analytic*csound*scal/dt`. Directly evaluating the quotient
derivative above loses accuracy when tiny w makes its two large terms
nearly cancel. The factored row retains every remaining divisor derivative;
the multiplier cannot relax the original energy residual tolerance.
The analytic local branch also uses `csound*scal/dt`. The
standalone checks also evaluate the ORIGINAL energy residual at convergence.
Positive-background LNA retains the original continuous equations; a
nonsingular row transformation at an exact equilibrium preserves the
eigenvalues. Do not use finite-dt predictors as physical LNA equations.

Boundary treatment:

1. Where the local residual has a common w factor under the proven conditions
   above, use its analytic divided expression and the dormant/positive
   `min(w, divided_expression)` row. Never divide by zero or add a denominator
   floor. With zero driving the selected row has a nonzero derivative.
2. When the full equation requires incoming energy, construct a positive
   trial w and keep active trial corrections positive. Extend
   `RSP2_adjust_vars_before_call_solver` to account for source, Eq, Lt and
   work/time-history terms. Predict new entropy moments jointly where they
   supply forcing. Accepted/start values remain unchanged; all physical
   coefficients and derivatives are current in the subsequent coupled solve.
3. Allow a zero row only when the original physical residual is zero there
   and the chosen dormant branch applies. Reevaluate with current neighbours
   and moment forcing. Finite accepted energy cannot be discarded to select
   dormancy. If a timestep has no admissible solution, a floor cannot repair it.
4. A fully quiet u-flag block needs its coupled viscous production considered.
   The check below has negative production diagonals but a positive collective
   mode and a positive solution of the original energy equations. A cell-local
   quadratic guess is not a complete activation rule. The predictor now sums
   all three Eq w derivatives together with the local source coefficient
   before its two directional sweeps. The earlier counterexample and two
   additional coupled blocks pass. This is not a general convergence proof.

`solver_support:Bdomain` now limits a decreasing active positive trial to at
least 0.1 times its current w. Forced cells and the analytic local branch can
reach exact zero. This bounds a Newton correction, not the physical w value.

Before this change, `solver_support:set_vars_for_solver` and
`hydro_vars:unpack_xh` replaced negative stored/trial w by
`RSP2_w_fix_if_neg` in the physical cache.
`auto_diff_support:wrap_w_00` still supplies a derivative of one. This can
make the value mapping and AD derivative inconsistent; the audit does not
establish how often this path occurs in a run. The proposed fix must keep
solver state and caches identical: repair admissible corrections before
evaluation, and reject invalid trials through the existing error path.
Invalid loaded/remeshed state needs an explicit error or justified state
reconstruction, not a cache-only floor. Do not change L_theta or ramp alfat.

Evidence in
[check_rsp2_shared_energy_boundary.py](check_rsp2_shared_energy_boundary.py)
and `output/review/rsp2_three_equation_20260919/shared_energy_boundary_checks.json`:

- 672 original-transport cases with unequal cells, zero/tiny receiving w,
  both face-weight options and theta_L=0.5/1: all converge within four updates;
  maximum mass-weighted energy error 3.14e-14 in dimensionless units.
- The narrow-cell example that stalls with the raw row and a square-root
  incoming-energy guess converges in four updates with the divided row.
- 162 scalar one-equation branch/root checks with signed production,
  quadratic storage/work/radiative coefficients and cubic dissipation.
- 720 two-cell three-equation cases, including zero/tiny w, variance-driven
  onset, u/v heating, alfat=0/0.01/0.2 and theta_L=0.5/1: all converge;
  maximum 16 iterations and original residual 8.99e-12. The original 84-case
  prototype still passes with this new treatment disabled.
- The collective u-flag onset counterexample identifies a remaining limit
  of the local predictor. Geometry/EOS/mean flow are fixed; signed
  time-centered viscous power, full hydro and MESA AD remain untested.

No MESA source edits, compilation, installation, evolution or LNA run were
performed for this shared-treatment investigation. This is a checked
candidate for positive cells plus an explicit boundary plan, not a claim
that all boundary cases or the full solver are now resolved.

### 2.5 Standalone verification of these decisions, 2026-09-19

Checked with [check_rsp2_three_equation_design.py](check_rsp2_three_equation_design.py),
using Python/mpmath/NumPy only, with results in
`output/review/rsp2_three_equation_20260919/design_checks.json`:

- 4,860 cases of the three temperature-row inversions versus 70-digit
  arithmetic; largest error scaled by the thermal-jump magnitude 4.87e-16.
- At a fixed floating-point neutral baseline, all three rearrangements retain
  a Y perturbation of 1e-20 that direct addition to gradL loses. This does not
  recover EOS digits already lost or prove arbitrary-resolution accuracy.
  That test prescribes a baseline making the cancellation exact; it does not
  derive the baseline from the actual pressure/temperature quadrature.
- Follow-up neutral-state check: 18 manufactured standard-row cases at
  three weights, three pressure jumps and Y=0 or 1e-8. Seventy-digit
  arithmetic confirms the mismatch described in section 2.1. At equal
  weights and pressure log jump 0.01, its bias is approximately 280 times
  the prescribed Y=1e-8 and reverses the inferred driving sign. Leading
  truncation terms are checked for both equal and unequal weights. This
  invalidates the previous entropy-gradient choice, despite its accurate
  inversion algebra. The replacement mapping checks are recorded below;
  actual MESA EOS/AD validation remains pending.
- Fixed-coefficient temperature-inverse derivative check: relative error
  2.23e-16. Full face-EOS/gravity/radiation-factor derivatives remain to be
  checked in the actual implementation.
- Local zero-state Jacobian has rank 2 with w, rank 3 with energy, for a
  nonresonant backward-Euler example. Existing Eq and Lt difference quotients
  nevertheless diverge in the energy coordinate at the specified boundaries.
- 1,000 random admissible-face cases satisfy the reconstructed cell covariance
  bound (largest ratio 0.908); the canceled local energy row reproduces its
  original residual with scaled difference below 1.57e-14.
- The empty-cell negative-source counterexample changes from -0.05 to zero
  under the chosen source reconstruction, without changing face Pi or Lc.
- Nine scalar imported-energy cases give positive initial guesses and positive
  roots of the full cubic residual for positive/zero/negative linear source
  coefficients. These are frozen-coefficient checks, not coupled u/v tests.
- The proposed scales are finite at zero turbulence, reproduce the chosen
  luminosity scale, and transform consistently with the entropy unit.

These checks verify design algebra and counterexamples. They do not verify
Fortran wiring, physical pulsation growth rates or stellar-model convergence.
No MESA compilation, installation, relinking or model run was performed.

The replacement mapper is checked separately by
[check_rsp2_three_equation_gradient_mapping.py](check_rsp2_three_equation_gradient_mapping.py),
with results in `output/review/rsp2_three_equation_20260919/gradient_mapping_checks.json`:

- 1,944 cases compare the mapped driving with the thermodynamic differential
  implied by the original temperature residual, using 80-digit arithmetic.
  These cover all three forms, dynamical on/off, signed composition offsets,
  unequal temperature/face weights, modified Hp, opacity floors, radiation
  factors, and distinct pressure/geometry time weights. Maximum relative
  difference: 7.48e-80.
- 1,296 cases preserve driving when gradL and Y shift at fixed gradT.
- 1,296 dynamical/actual-log neutral cases give exactly zero in double
  precision. 3,888 signed small-Y cases down to 1e-20 agree with the
  high-precision result to 2.94e-15 relative error.
- 216 complex-step derivatives, including variable face/opacity/radiation
  coefficients, agree with differentiation of the uncanceled differential
  to 1.71e-15 when scaled by max(1,abs(reference derivative)). These are
  Python checks; they do not exercise Fortran AD.
- 18 static-HSE cases preserve tiny driving in high precision. Rounding the
  cell pressures leaves up to 1.58e-14 of the neutral pressure-force term
  in double precision. The static path deliberately retains this mismatch.
- The actual-log expression remains finite for zero and reversed pressure
  gradients; it never divides by that gradient.

The tests manufacture thermodynamic coefficients; they do not call the EOS,
exercise real mesh boundaries/reconstruction, or validate stellar convergence.
`constant_L` retains a measured differential and its cell-value precision
limit; it does not share the analytic neutral-cancellation guarantee.

A second standalone check,
[check_rsp2_three_equation_zero_boundary.py](check_rsp2_three_equation_zero_boundary.py),
solves a two-cell, four-variable system (two w values, one face Pi and one
face entropy variance) on unequal masses. It includes the actual cell/face
velocity dependence of the two viscous-heating placements, while holding
their geometric/mean-flow coefficients fixed. Results:

- 84/84 selected cases converged with nonnegative w and entropy variance.
- Cases include alpha_t = 0, 0.01, 0.2; alpha_m = 0, 0.25; positive/negative
  entropy driving; dormant cells; zero initial flux; and finite variance
  with initially zero kinetic energy.
- Largest ORIGINAL energy/moment residual at convergence: 8.33e-12 in
  the prototype's dimensionless units. Both branch and original residuals
  are checked; projection of a trial correction alone cannot pass a case.
- The prototype needed at most six iterations for these cases. This says
  nothing about iteration counts for the full MESA equations.
- It does not test evolving velocities/time-centered strain reversal,
  momentum-work cancellation, reconstructed EOS derivatives, remeshing,
  the full spatial mesh, MESA line searches, or LNA.

Full case records: `output/review/rsp2_three_equation_20260919/zero_boundary_checks.json`.
These branch decisions are now wired into the source. No physical closure
or covariance bound guarantee is inferred from the standalone checks.

## 3. Complete infrastructure map

| Area | Existing source/routine | Required change and invariant |
| --- | --- | --- |
| Controls | `star/defaults/controls_dev.defaults`, `star_data/private/star_controls_dev.inc`, `star/private/ctrls_io.f90` | Default, declaration, namelist, copy-in and copy-out for one switch and two coefficients. |
| Active state | `star_data/public/star_data_step_input.inc`, `star/private/init.f90` | A stored `RSP2_3equation_flag` distinguishes the loaded layout from the requested control; initialize false. |
| Indices | `star/private/alloc.f90:set_var_info` | Add two variables after Y only in active RSP2 three-equation mode; two named residual rows at the same indices. |
| Work arrays | `star_data/public/star_data_step_work.inc`, `alloc:star_info_arrays` | Two physical moment caches and two fixed start-state scale arrays; allocate/copy/free/size-check through existing do1 machinery. |
| AD slots | `star/private/auto_diff_support.f90`, `star_utils:unpack_residual_partials` | Use spare triples i_xtra1/i_xtra2; named wrappers; map all m1/00/p1 derivatives to new columns. Keep AD dimension unchanged. |
| State unpack | `hydro_vars:unpack_xh` | Load caches from authoritative xh only when their indices exist. |
| Trial updates | `solver_support:set_solver_vars` | Populate caches from xh_start+solver_dx, including numerical partial tests; xh alone is not current during trial evaluation. |
| Accepted state | `struct_burn_mix` initial-guess copy and `hydro_solver_step` | Copy new caches into xh for guesses; successful solve already copies all active columns. |
| Residual dispatch | `hydro_eqns:eval_equ_for_solver` | Call both new rows under the RSP2 branch, with the same nvar guards/error handling as adjacent rows. |
| Moment physics | `hydro_rsp2` | Shared rhs routine for nonlinear solve and LNA; no duplicate LNA closure. |
| Lc | `hydro_rsp2:compute_Lc_terms` | New mode uses Pi directly. False branch retains old w*PII expression. |
| Source/Dr | `compute_Source`, `compute_Dr_div_w`, `compute_C` | Covariance-based source; no direct Dr in new mode; same C in all energy consumers. |
| Zero-w branch | `do1_turbulent_energy_eqn`, `RSP2_adjust_vars_before_call_solver` | Use the reconstructed source, canceled local branch where factorization is valid, full equation for imported Lt/Eq, and explicit dormant constraints; see section 2.4. |
| Initialization | `set_flags:set_RSP2_flag` and new mode-state helper | Insert/remove only the two moment columns when RSP2 was already on; do not re-seed w or Y. |
| Load/startup | `star/job/run_star_support.f90`, `read_model:finish_load_model` | Reconcile requested mode after loaded state/EOS are valid and before startup LNA; avoid module cycles. |
| Step transition | `evolve:prepare_for_new_step` | Reconcile control at a step boundary before mesh snapshots/new_generation; never during Newton or retries. |
| Photo state | `star/private/photo_in.f90`, `photo_out.f90`, `star_data_def.inc` | Preserve stored layout before reading xh; keep version-20 compatibility through a versioned extension. |
| Saved models | `read_model.f90`, `write_model.f90` | Dedicated file_type bit and two columns; preserve old model decoding and column counts. |
| Normal remesh | `mesh_adjust:do_RSP2_face_var`, `do_etrb` | Reuse shape-preserving face interpolation for new moments; retain conservative cell e_t remap. |
| Split/merge | `adjust_mesh_split_merge.f90` | Generic shifts already carry xh; initialize newly created face values, keep surviving faces on merge, refresh caches. |
| Envelope remesh | `tdc_hydro_support:remesh_for_TDC_pulsations` | Use existing face interpolation for both new slots before resizing/repacking; preserve w^2 overlap remap. |
| Profile output | `star/private/star_profile_def.f90`, `profile_getval.f90`, `star/defaults/profile_columns.list` | Register IDs/names and values for physical moments and covariance diagnostic; return defined values in one-equation mode. |
| Mixing diagnostics | `star/private/mix_info.f90`, `turb_info.f90` | Audit existing w/Lc/Lt consumers and mixing classification; do not introduce a second mixing law from Pi. |
| LNA map | `star_LNA_support:setup_star_LNA_var_map`, `equation_id_for_star_LNA` | Two additional variables and equations only for active new mode. |
| LNA rows | `assemble_rsp2_turbulent_rows`, `add_ad_partials_to_A/B` | AD of the shared rhs goes to A; identity moment inertia to B; forced faces are algebraic zero rows. |
| LNA partition | `star_LNA_support:star_LNA_var_is_dynamic` and reduced matrix assembly | Register both moments as dynamic candidates; classify boundary rows by their actual B rows. Adding matrix columns alone is insufficient. |
| LNA luminosity | `star_LNA_turbulence_closures:rsp2_convective_luminosity_for_star_LNA` | Use covariance directly; keep existing frozen/perturbed flux policy. |
| LNA reports | eigenfunctions, matrix row names, setup diagnostics | Name/output the two perturbations and include moment equilibrium defects; no silent unused columns. |

## 4. Activation and state lifecycle

Control and active flag must not be used interchangeably. The file's active
flag determines how many xh rows to allocate before reading its data.
The control is applied only after that read succeeds.

Activation sequence for an existing one-equation RSP2 model:

1. Finish reading the original layout and compute its current EOS/Y/Lc.
2. Save old indices/counts; insert two hydro columns after Y using MESA's native
   set_flags/update_nvar_allocs pattern. Preserve other hydro columns,
   including rotation variables; abundances are in the separate xa array.
3. Initialize Pi from the old physical Lc, then the missing variance.
4. Initialize new xh_start values and invalidate unavailable previous-moment
   history rather than pretending zero historical moments are physical.
5. Refresh caches and derived variables. Preserve w, L and physical gradT.
   Convert Y to `gradT_old-gradL_new` when the neutral reference changes.
6. Record the active flag; subsequent set_vars and restarts preserve moments.

For initial activation of RSP2, its existing w/Y initialization and optional
remesh happen normally, then the same new-moment conversion applies.
For disabling RSP2, remove the moment columns before removing Y and w.
For toggling only the three-equation control off, remove just its two columns.
Repeated calls with matching active/requested state must do nothing.

Implementation should not require another public star_job control. Prefer
extending the existing star_set_RSP2_flag entry to synchronize its selected
mode, with a dedicated internal helper if this keeps the existing w/Y
transition intact. Call synchronization from startup and a new-step boundary.

### 4.1 Photo/model format strategy

Photo version 21 adds an active mode record after the original initial
flag record. The reader accepts both 20 and 21;
version 20 implies no stored moments. The old indices record can remain
unchanged: new moment indices are derived by set_var_info from the stored
mode flag. xh's generic record then carries the new values exactly.

Saved models use file_type bit 17 to identify the two extra moment
columns. Verify read1_model's nvar/extra-column accounting rather than adding
two extra counts twice. Preserve old one-equation .mod files. On load,
original .mod/photo layout is authoritative even if the new control differs.

Restart the same mode: preserve Pi/variance exactly. Change mode: initialize
only missing moments once. Retry: restore xh and regenerate caches through
unpack; never run the mode initializer during retry/redo.

Verified details from the actual lifecycle:

- `alloc:update_nvar_allocs` resizes xh/xh_old/xh_start, solver arrays and
  prev_mesh arrays and sets `prev_mesh_species_or_nvar_hydro_changed=.true.`.
  It does not perform the physical hydro-column permutation for the caller.
- `evolve_support:new_generation` copies every active hydro row into xh_old;
  `set_current_to_old` restores every row. No moment-specific history arrays
  are needed if the active layout stays fixed throughout an attempted step.
- `evolve:prepare_to_retry` can first restore `prev_mesh_xh` into xh_old.
  Synchronize layout before taking that snapshot; never restore an older
  layout after activating the new variables. Use nz_old when touching old
  arrays, rather than assuming nz_old equals nz.
- Startup runs `do_star_job_controls_after`, optional envelope remesh,
  `extras_startup`, then `maybe_do_star_LNA`. The synchronization must precede
  the first consumer requiring the selected equations, including startup LNA.
- Saved model bit 17 now identifies the moment layout.
  `read1_model` sizes its vector using nvar_hydro plus existing file-only
  extras; adding two solver variables already enlarges that capacity by two.
  Check column consumption with nvec, including optional previous-model data;
  do not add an extra +2 on top of the enlarged hydro count.

### 4.2 Startup seed and its honest limitations

When creating missing moments at an unforced, unstable face with positive
old convective luminosity, use

```math
Pi_f = Lc_old/[4 pi r_f^2 rho_face,new T_face,new],
<s'^2>_f = 3 Pi_f^2/[2(e_t)_f],  (e_t)_f>0.
```

Only when the old and new face conversion coefficients agree does this equal
`average(w)_old*PII_old`. Stable and neutral faces, faces with nonpositive
old flux, and forced nonturbulent faces use zero moments. Retain cell
turbulent energy in every case. A retained nonzero flux with zero face
energy is inconsistent and is rejected without a division floor.
Compute the target gradL and unpack Y=gradT_saved-gradL before classifying
the face. This preserves the imported heat flux only at the selected
unstable faces and adopts a fully correlated MLT parcel variance there.
It does not reconstruct unknown time history or guarantee zero new
residuals. Existing three equation restart states bypass this initializer.

A quiet, static local model has the analytic equilibrium derived in the
proposal; its no-radiation limit reproduces the one-equation local flux law.
Do not silently replace an evolving model's w by that equilibrium. A wholly
zero local three-moment state is an absorbing state without imported
energy/perturbations; the local canceled row alone is not a physical onset source.
Any future onset prescription must be explicit, not hidden in an initializer.

## 5. Energy/work audit matrix

Verified `hydro_energy:get1_energy_eqn` and its
`setup_d_turbulent_energy_dt`: BOTH dedt and eps_grav residuals already
subtract d(e_t)/dt. The eps_grav row also uses the simple work path with
thermal pressure work excluded because eps_grav already contains it.
Neither row needs a new explicit -C term. The exchange is recovered by
combining it with the turbulent-energy row; adding -C separately would
double-count it. The new physics enters through C in the turbulence row and
Lc in the face flux balance, which determines L in the energy divergence.

`hydro_energy:unpack_res18` calls the shared `unpack_residual_partials`, despite
its historical name. The new AD columns therefore belong in that shared
mapping, not in a separate rewrite of the energy Jacobian.

| Velocity | Energy form | Work choice | Required invariant |
| --- | --- | --- | --- |
| v_flag | dedt | simple PdV | New moment exchange does not create or destroy total energy. |
| v_flag | dedt | conservative work | Existing Eq/Uq pairing retained; no new wider derivative loss. |
| v_flag | eps_grav | both work settings | Preserve d(e_t)/dt and the existing forced simple-work path; no added -C. |
| u_flag | dedt | both work settings | Existing cell momentum/viscous power retained. |
| u_flag | eps_grav | both work settings | Same accounting; RSP2-supported eps_grav path remains active. |

Check nonzero alfat separately: Lt appears consistently in total luminosity
and turbulent-energy divergence. The two local entropy moment closures do
not disable kinetic-energy transport or its conservative remap.

## 6. Mesh and boundary invariants

- Face moments are interpolated as point values; do not overlap-conserve
  their dq-weighted sums as though they were cell turbulent energy.
- Use existing monotonic piecewise interpolation for normal/envelope remesh.
  Signed Pi remains signed; positive variance must not acquire interpolation
  undershoots. Preserve unchanged faces exactly where the current Y routine
  already does so.
- A split creates one interior face: interpolate from the original bounding
  faces with the child mass fraction. Both old boundary faces survive.
- A merge removes an interior face: surviving-face moments are retained.
  Do not average them into a new cell-centered moment.
- Generic xh shifts do not automatically update physical caches. Verify
  that remesh exits run unpack/set_vars before any moment physics consumes
  them, or explicitly update caches where adjacent Y/w code does so.
- Forced surface/inner nonturbulent faces need algebraic Pi=variance=0 rows
  and matching zero heat flux. Clarify excised boundaries versus full center.
- Monitor the covariance bound after remap, but do not silently clip Pi to
  impose a new flux limiter. Interpolating different moments separately does
  not prove preservation of their joint bound.

## 7. LNA integration

Use one nonlinear rhs implementation for both moment rows. With exp(sigma*t),

```math
sigma delta Pi = delta rhs_Pi,
sigma delta <s'^2> = delta rhs_s2.
```

The matrix map gains two variables per active zone. All six spare AD indices
must map to their correct zone/variable. Source and luminosity perturbations
must include both new moments where appropriate; old dPII/dY is no longer
an independent closure for this mode. Keep radiation and Lt perturbations,
work decomposition, boundary conditions and u/v distinctions intact.

The current setup requires `star_LNA_perturb_turbulent_energy=.true.` for
RSP2; retain that requirement for all three moment rows. The independent
`star_LNA_perturb_convective_flux` option still controls delta Lc in the
energy/flux closure. Do not infer that this switch removes the covariance
evolution equations. Preserve its existing frozen-flux meaning explicitly.

For LNA's inner cut, deeper background values remain fixed according to
existing map rules. Do not force an artificial zero-moment condition at the
cut merely because map%nz differs from the physical model's nz.

LNA diagnostics must report static moment residuals: a switched pulsating
snapshot is not automatically a suitable equilibrium background. No change
to the published current-model PDF until the new implementation is checked;
then distinguish optional new-mode equations from established one-equation
RSP2. Do not describe proposed behavior as already validated.

## 8. Verification plan without compiling or running MESA

1. Inspect every control occurrence: default/declaration/namelist/two copies.
2. Verify active-state initialization, flags, index counts/names and AD
   column mapping; aliases must not collide with existing active derivatives.
3. Parse changed Fortran and run Fortitude; these are not compilation.
4. Standalone numerical/complex-step checks of both rhs equations, buoyancy
   source and face luminosity at unequal zones, signed flux and weak energy.
5. Check the local stationary solution, including radiative timescale factors
   one and two; demonstrate original local closure recovery at default alfa coefficients.
6. Check gas/turbulent exchange cancellation and telescoping luminosity
   divergence with nonuniform dm; preserve existing Eq/Uq checks.
7. Check activation/removal hydro-column permutations for u/v/rotation
   variants and differing old/current nz; preserve the separate xa array and
   the hydro/chemistry offset used in the combined solver vector.
8. Check .mod/photo layout symmetry and version-20 read path statically and
   with standalone record fixtures where practical, without a MESA executable.
9. Check normal, split/merge and envelope interpolation using unequal zones
   and signed/positive data, including unchanged-face and boundary cases.
10. Check LNA derivatives/inertia/AD mapping against the nonlinear moment
    equations; count rows and variables and check forced faces.
11. Run existing `notes/check_rsp2_independent_Y.py` and
    `notes/check_tdc_cell_mixing_length.py`; distinguish intended assertions
    needing extension from actual regressions.
12. `git diff --check`, inspect final scoped diff, update this plan with actual
    outcomes and remaining runtime tests. No compile/install/run unless the
    user explicitly authorizes that next action.

Runtime checks for the user later: old photo with flag off; old photo with
flag on; new photo restart; new .mod load; both energy forms and velocities;
alfat zero/nonzero; each temperature-gradient form; normal remesh and AMR;
short time-step convergence study and LNA comparison on a static background.
Passing source/algebra checks is not a claim those runtime checks passed.

## 9. Implemented source and checks

| Area | Current state |
| --- | --- |
| Controls and names | `RSP2_alfa_pi` and `RSP2_alfa_phi` default to one; source constants carry the local MLT calibration. |
| Active state and layout | Two columns after Y; transition before startup LNA and each new step; unchanged mode returns without conversion. |
| Moment rows | Both residuals dispatched, common rhs used by LNA, fixed start scales and complete AD maps. |
| Heat and energy exchange | Pi supplies Lc and reconstructed buoyancy. Dr is absent only in three equation mode. Existing total energy and work terms retained. |
| Zero turbulent energy | Analytic local branch, positive full energy row, coupled heating/import guess, consistent trial domain and no cache substitution. |
| Saved state | Version 21 photos, version 20 reader, model bit 17, generic hydro history and retry copying. |
| Mesh | Monotonic face interpolation for normal/envelope remesh; point values on AMR; conservative cell energy remap retained. |
| LNA and profiles | Two moment variables and dynamic rows, dormant constraints, instantaneous gradients, eigenfunctions, equilibrium defects and three profile columns. |
| Validation | Standalone numerical/source checks and Fortitude pass. Binary restart, actual MESA interpolation, Fortran AD, full hydro and eigenmodes remain untested. |

### Audit evidence recorded during the planning pause

- No additional Fortran/default/include files were edited during the pause.
  At that point the provisional diff was 13 tracked files, 280 insertions,
  7 deletions. The subsequent gradient-mapping edits are recorded below.
- Read the actual new-generation, retry, redo and pre-remesh snapshot paths;
  verified generic hydro-row copying and the layout-change invalidation flag.
- Read all three energy residual branches; verified d(e_t)/dt is present in
  each and that energy AD unpacking uses the shared mapping.
- Located all three profile registration files and the explicit LNA dynamic
  variable classifier, in addition to its variable/equation/AD maps.
- The spare AD triples are slots 28--30 and 31--33. Repository search found
  only the generic zero-valued placeholder wrappers outside the provisional
  implementation; w_div_wc uses slots 22--24, so those are not the same slots.
- The pause audit found two trailing spaces in provisional `alloc.f90`
  additions. Both were removed when source editing resumed.
- MESA compilation, installation, executable runs, and case relinking remain
  unperformed. No numerical-convergence claim follows from this source audit.

### Source work resumed for the gradient mapping

The new support module is listed in `star/Makefile`. QHSE/gravity and dPrad
factor calculations are shared to avoid duplicating the temperature-row
coefficients; existing callers retain their default time centering and
radiation factors. The optional instantaneous argument prepares the helper
interface for LNA; it does not implement the new LNA rows.

The mapper sets gradL and thermal driving together, uses the existing
composition term explicitly, and is selected only by the inactive
`RSP2_3equation_flag`. That flag is still only initialized false. No control,
model/photo conversion or equation dispatch has been added in this step.

Source review caught and corrected bad-number checks that initially passed
an AD value to `utils_lib:is_bad`; they now pass `%val` like existing MESA
callers. The new module passes Fortitude and `git diff --check` is clean.
The shared-module dependency check finds no cycle reachable from the
temperature, momentum, RSP2 or new support modules. No compilation or model
execution was used to obtain these results.

## 10. Execution checklist

- [x] Preserve unrelated source, notes, cases and outputs.
- [x] Commit the shared gradient changes separately.
- [x] Implement row consistent entropy driving and check its algebra and derivatives.
- [x] Wire common face coefficients, cell source reconstruction and start scales.
- [x] Replace shared zero w residual/domain handling and check coupled heating examples.
- [x] Wire controls, mode conversion, photo and model state.
- [x] Wire normal, envelope and split/merge mesh paths.
- [x] Wire LNA variables, rows, closures, inertia, diagnostics and profiles.
- [x] Check source dependencies and parse changed Fortran with Fortitude.
- [x] Run standalone energy, gradient, moment, viscosity and remap regressions.
- [x] Finish manuscript rendering and final diff review.
- [x] Compile and install MESA after the user's explicit authorization.
- [ ] User checks actual Fortran AD and model evolution.
- [ ] User verifies photo/model round trips, retries and each remesh path.
- [ ] User verifies LNA on a relaxed static background and timestep convergence.

The remaining numerical validation requires MESA execution, which only the
user is authorized to perform in this task. The test matrix in section 8
specifies the cases. Source checks do not replace those runs.

## Implementation progress after the gradient commit

The gradient helper commit is `1700f15ac`. Its message follows the requested
MESA/EbF wording. Subsequent source changes remain uncommitted while the full
mode is being completed and checked. No MESA compilation or model run has
been performed.

The user selected `Pi` for `<v_r*s>` and `Phi` for the full `<s*s>` variance.
The stored fields, indices, controls and AD wrappers now use these names.
`PII` retains its original meaning. The MESA implementation skill records
this naming preference and concise MESA/EbF comments, notes and commit text.

Source integration completed in this continuation:

1. Common face thermodynamics in Lc and buoyancy, the resolved pressure
   difference, and the cell w factor in the reconstructed source.
2. Pi and Phi residual dispatch and luminosity based row and column scales.
3. Mode conversion preserving the old Lc and physical gradT, with unavailable
   moment history invalidated. Synchronization occurs after startup loading
   and before a new step, never during retry or Newton iteration.
4. Photo version 21 with version 20 reading retained. Saved models use bit 17.
5. Normal, envelope and split/merge remeshing for the two face variables.
   The existing conservative remap of cell turbulent energy is retained.
6. LNA variable and AD maps, moment inertia and shared continuous rhs, plus
   the instantaneous gradient mapping in the temperature and luminosity rows.
7. Pi, Phi and covariance excess profile columns, moment eigenfunctions and
   dimensional moment equilibrium residual reporting.

The v velocity strain uses the face interpolation of the two neighboring
cell slopes `(v_outer-v_inner)/(r_outer-r_inner)`. This gives the exact slope
of an affine radial velocity on an unequal radial grid. For cell velocity u,
the face strain uses the existing mass to radius conversion times the cell
velocity difference. Both use current coefficients and AD derivatives.

The shared energy row now retains its AD divisor at positive w. The local
zero start branch is canceled analytically; forced cells can reach zero.
The initial guess includes source, damping, Eq, accepted/current Lt and
pressure work, with two directional sweeps. A quiet u region uses the sum
of neighboring Eq derivatives as well as local production. Remaining Phi
at zero kinetic energy seeds a trial Pi/w response without changing accepted
state. These are initial guesses, not changes to the physical equations.

Negative stored or trial w is rejected. Pi is signed; Phi must be finite and
nonnegative. The old `RSP2_w_fix_if_neg` control is accepted for compatibility
but no longer substitutes a cached value. No luminosity time weights or
physical Eq/Uq terms were changed in this implementation.

Latest standalone evidence, in
`output/review/rsp2_three_equation_20260919/implementation_checks.json`:

- 240 moment/source derivative cases, including a zero energy cell beside a
  finite energy face; maximum scaled difference 4.95e-13.
- 100 continuous/discrete moment block comparisons.
- Three quiet coupled heating blocks, including the earlier counterexample;
  largest original energy residual 2.55e-21.
- 160 transport cases using the new directional guess; all converged in at
  most four Newton updates in the isolated fixed coefficient problem.
- 500 AMR conservation/point interpolation and affine velocity cases.
- Control, AD, row, column and file layout source checks; no binary round trip.
- No module cycle in the selected star source graph.

The earlier shared energy suite also passes 672 transport, 162 local branch
and 720 coupled two cell moment cases. The gradient mapper, independent Y,
cell mixing length/viscosity and conservative remap suites pass. The renamed
face remesh helper required an expected name update in the source regression;
its interpolation and conservation checks were retained.

The implementation adds no realizability limiter. Full hydro convergence,
physical calibration, actual mesh interpolation/Fortran AD and eigenvalues
still require the user tests. No compilation or installation has occurred.

Final source review passes `git diff --check` and Fortitude on all 25 changed
Fortran files. The manuscript was rebuilt and the added equation and closure
pages were visually checked. It has 35 pages, no overfull boxes or unresolved
references, and no text outside the checked page margins. The source keeps
the new coefficient controls as `RSP2_alfa_pi` and `RSP2_alfa_phi`, both one;
their fixed factors and the alfar/alfat normalization are documented. Further
style review is deferred at the user's request. The mode remains off by
default and the implementation changes remain uncommitted.

## Installation and separate RSP3 manuscript, 2026-09-19

The user explicitly requested installation before the manuscript so they can
relink and test. Initialized MESASDK_ROOT=/Applications/mesasdk, then set
MESA_DIR=/Users/owner/Documents/Software/dev/test/mesa, unset GYRE_DIR and
ran ./install with NPROCS=10. Installation and standard package checks passed
with exit status 0. The hydro_rsp2, hydro_gradient_support, star_LNA and
set_flags objects match their installed libstar.a members byte for byte.
Logs and hashes are in output/build/rsp3_install_20260919.log and .json.
No user case was relinked or run. This establishes compilation and package
checks, not convergence or physical calibration of the new mode.

The requested separate manuscript is notes/rsp3_manuscript.tex and .pdf.
It formalizes the equations, residuals, linearization, normalization and
state/mesh handling. Development progress stays here; implementation status
text was removed from the LNA manuscript. Its scientific equations remain.

Coefficient clarification, corrected after the startup audit: the unit
Pi/Phi controls multiply the local MLT constants 6*sqrt(2/3) and 4*sqrt(2/3).
Braun et al. Table 1 instead uses alpha_Pi=2.155 and alpha_Phi=2.0 for
nonlocal spatial transport. Dividing these numbers by the local constants
does not give calibrated RSP2 decay controls; the earlier mapping was
misleading and is withdrawn. Alfar may remain zero,
which removes radiative cooling of the two moments but retains turnover
relaxation. No defaults were changed in response to these questions.

The RSP3 manuscript is 12 pages, compared with 35 for the LNA manuscript.
It pairs the three turbulence equations with their finite step residuals
and static linearizations, and retains both energy forms, both velocity
grids, temperature row mapping, Eq/Uq work, source reconstruction, zero w
handling, initialization, remeshing, saved state and LNA. The manuscript
does not carry development status. The LNA manuscript was rebuilt after
removing its implementation status sentences; its scientific content was
retained.

All 12 RSP3 pages and the affected LNA pages were rendered and visually
reviewed. Both PDFs have no overfull boxes, unresolved references or text
outside the checked page margins. PDF checks are recorded in
output/review/rsp3_manuscript_20260919/pdf_checks.json.
The document review checked energy residual signs against hydro_energy
and recorded the existing conservative energy Jacobian stencil limitation.
No Fortran source or user case controls were changed during this manuscript
and installation step. The user still handles relinking and model runs.

## Structural timestep correction after the first user run, 2026-09-19

The user's attached trace and current case history show 83 accepted models,
zero retries, and termination at `min_timestep_limit`. After nine Newton
iterations in model 1, the solve takes two or three iterations. Nevertheless,
every accepted timestep is reduced by exactly 0.8:

```math
dt_n = 100 (0.8)^(n-1) seconds.
```

The last accepted dt is 1.130782121458171e-6 seconds and the requested next
dt is 9.046256971665369e-7 seconds, below the 1e-6 second minimum. The run
advances only about 500 seconds. This is a timestep controller failure
after convergence, not repeated Newton failure.

The integration omitted Pi and Phi from the exclusions in
`timestep:eval_varcontrol`. The existing structural change norm deliberately
excludes luminosity, velocity, RSP2 w and Y, angular momentum and other
nonstructural variables. It was instead adding smoothed raw changes of the
new dimensional moments to changes of logarithmic structure variables:

```math
varcontrol = sum_j sum_k |smooth(xh_new(j)-xh_old(j))|
             / (nterms * sum_j max(1,|xh_old(j,1)|)).
```

The usual density factor and endpoint weights are omitted from this display.
Both surface moments are zero, so each erroneous moment term adds only one
to the denominator scale despite very large interior cgs amplitudes.
The first profile has Pi of order 1e15, reconstructed from the diagnostic
PII times its face w. Its raw dimensional change has no place in that
structural norm. In `do_timestep_limits`, the `max decrease` label is used
when this structural limiter requests a decrease beyond
`min_timestep_factor`, which defaults to 0.8.

The source correction adds `i_Pi` and `i_Phi` to that existing exclusion
list. It changes no residuals, initialization, time weights, controls,
energy conservation or moment solver tolerances. It retains the ordinary
thermal and hydrodynamic timestep limits. A dedicated moment accuracy
criterion, if wanted later, would require a separate dimensionless measure;
raw cgs moments are not one.

The standalone regression in
`check_rsp2_three_equation_implementation.py:check_timestep_varcontrol`
first failed against the uncorrected source. It compares identical thermal
changes with and without moments, both velocity layouts, optional trailing
variables and rescaled moment units. The corrected source must give the
same structural measure and zero measure for moment-only changes.
The corrected source passes all 18 targeted cases, Fortitude and
`git diff --check`. Reinstallation follows the existing install authorization,
with MESA_DIR explicitly set to this checkout. The user retains case relinking
and model execution. Installation outcome is recorded below.

Installation completed successfully with exit status 0 and
`MESA_DIR=/Users/owner/Documents/Software/dev/test/mesa`. The package checks
passed, and the installed `libstar.a` member `timestep.o` matches the newly
compiled object byte for byte. The build log and hash verification are in
`output/review/rsp3_timestep_20260919/install.log` and `install_check.json`.
No user case was relinked or run. Repeating the startup run remains the
runtime validation of this correction.


## Authorized 2000 step run, 2026-09-19

The user authorized running `dev_TDC_RSP2_Cepheid` for 2000 steps. The
existing executable is newer than the installed library. Run it directly
with this checkout as MESA_DIR and eight threads. Set max_model_number=2000
and profile_interval=1; add Pi, Phi and Pi_covariance_excess to the profile
columns. Keep all physics, initialization, solver and timestep controls.
Preserve previous outputs and exact configuration in
`output/review/rsp3_2000_steps_20260919/before`. Restore the two edited
configuration files after completion. New run products stay in the normal
case output paths. The run manifest records the executable hash and changes.

- [x] Preserve existing outputs and configuration.
- [x] Complete the run or record its actual early termination.
- [x] Measure accepted timestep and iteration history, retries and limiting rows.
- [x] Trace persistent problems through moment, energy and gradient equations.
- [x] Restore temporary output and stopping controls; record conclusions.

The run completed 2000 accepted steps with 427 retries and 2000 profiles.
Models 21 through 2000 advanced only 36.7077 seconds. Accepted timesteps
remained between 0.0134435 and 0.0271656 seconds after model 20. 1981 steps
took ten iterations. The ordinary Phi residual dominated 30438 of 30520
printed iterations; it is not the previously corrected varcontrol problem.

At face 145, Phi reaches zero while Pi stays positive and -ds/dr is negative.
The variance equation then demands dPhi/dt=2*(-ds/dr)*Pi<0, but Bdomain clips
the negative Phi update. Reconstructing the final row from profiles gives
9.54099e-5, matching the printed 9.541e-5. It passes only after the residual
tolerance relaxes at iteration ten. Timestep growth repeats the failure.

The imported model has w_start=7.59294e-8 cm/s at this cell. The correlated
seed Phi=1.5*Pi^2/etrb_face can nevertheless create a large entropy variance
because the inherited flux is proportional to w; that factor cancels in
the seed. Review this initialization and the RSP-to-RSP2 gradient conversion
ordering before deriving imported moments from the recomputed old flux.
Separately, a homogeneous continuous ODE check demonstrates that the local
moment closure can cross Phi=0 with dPhi/dt<0 in stable stratification, with
both the user's coefficients and unit defaults. This remains a physical
closure problem even without a mesh or Newton solver. A safe correction
is not established by these diagnostics. Do not claim runtime readiness.

The full equations, source references, statistics, and reproducible analysis
are in `output/review/rsp3_2000_steps_20260919/analysis.md`, `analysis.json`,
and `analyze.py`. The plot is `run_summary.png`. All temporary configuration
changes were restored byte for byte; the original outputs are preserved in
that review directory's `before` subdirectory. New outputs remain in the
case. No source or physics settings were changed during this test.

## Seed and local coefficient analysis, 2026-09-19

The user asked whether Phi needs a seed. A separate positive Phi floor is
not required when the face kinetic energy is nonzero. With Pi=Phi=0,

```math
dPi/dt = (2/3)*etrb_face*(-ds/dr),
dPhi/dt = 0,
d2Phi/dt2 = (4/3)*etrb_face*(ds/dr)^2.
```

Thus the resolved entropy gradient generates flux and then positive
variance. Both signs of the entropy gradient give a positive leading
variance. These are initial continuous derivatives, not a replacement
for the implicit residuals. If all three moments are exactly zero there
is no spontaneous local onset; kinetic seeding or imported energy remains
a separate issue.

There is an additional coefficient incompatibility in the current run.
For a homogeneous stationary convective layer without radiation, shear
or spatial transport, eliminating Pi and Phi from the three source
equations gives

```math
Pi^2 / [(2/3)*etrb_face*Phi]
  = (RSP2_alfa_phi + 2*RSP2_alfad)/(3*RSP2_alfa_pi).
```

Unit controls give one. The user's controls (alfa_pi,alfa_phi,alfad)
=(0.4,0.6,1) give 13/6. Even an exact local stationary initializer with
these controls exceeds the covariance bound under the adopted isotropy
assumption. This condition concerns the stated local, nonradiating limit;
it is not a universal stability criterion or a statement about all of
Braun's nonlocal solar models. Her fiducial transport coefficients cannot
be presumed calibrated local decay rates for this pulsation model.

The standalone script
`output/review/rsp3_2000_steps_20260919/check_moment_seeds.py` verifies both
stationary solutions and integrates quiet stable-layer seeds with
Pi=Phi=0. Four initial energy scales and two tolerances were checked.
For unit controls all eight trajectories remain nonnegative over the
50-unit test interval. For controls 0.4 and 0.6 all eight reach the
variance boundary with a nonzero flux. The output is
`moment_seed_checks.json`. These are dimensionless homogeneous tests,
not additional MESA runs or a proof of general positivity. The earlier
counterexample with entropy variance initially present still applies
even with unit controls.

The first controlled startup candidate is therefore: correct target
gradL before converting Y; initialize the missing moments in the nearly
quiet stable layers without manufacturing a finite entropy variance
from the inherited velocity floor; and establish a baseline with the
local unit coefficients. With retained positive w, Pi=Phi=0 is a physical
choice for those quiet layers and the equations generate subsequent
fluctuations. This initialization changes the initially assumed heat
flux and must not erase actual moment history, countergradient flux in
an evolved three equation model, or restart state. Convective layers
require a compatible coupled local equilibrium or a documented physical
transient assumption, with luminosity and temperature consistency checked.

A successful startup comparison would establish improvement for this
case, not remove the general transient covariance issue. No source or
case controls were changed for this analysis.


## Startup correction implementation, 2026-09-19

The user authorized the proposed startup fixes. The invariant is to retain
cell turbulent energy and the physical temperature gradient while creating
missing face moments. Stable faces start with Pi=Phi=0; an unstable face
with positive imported Lc retains that flux and its correlated variance.
This is an explicit startup assumption, not a projection during evolution.
Same-mode restart and retry state remain untouched.

- [x] In set_RSP2_flag, save gradT before changing layouts, compute target
      gradL with Y=0, assign Y=gradT_saved-gradL, then refresh the old flux.
- [x] In set_RSP2_3equation_flag, update and unpack Y before initializing
      moments so get_rsp2_thermal_gradient sees the preserved gradient.
- [x] Gate init_rsp2_moments on the actual entropy driving, without a
      velocity floor; retain the existing energy and forced-face rules.
- [x] Set only the two new coefficient controls to 1d0 in the test case.
      Clarify local decay versus nonlocal transport coefficients in docs.
- [x] Check stable/unstable initialization, gradient conversion, restart
      guard, source wiring and compilation. Update the rsp3 manuscript.

The before-state files are saved under
`output/review/rsp3_startup_fix_20260919/before`. This work does not claim
to enforce realizability of the continuous moment closure.

Validation is recorded in
`output/review/rsp3_startup_fix_20260919/standalone_checks.json`:
27 gradient conversions and 45 moment seeds pass, including the same-mode
restart guard and initialization with refreshed Y. The full existing
standalone suite also passes. These are algebra and source checks, not
executed MESA restart or coupled hydro tests. Fortitude passes for both
changed Fortran files, and `git diff --check` passes.

Installation completed successfully with
`MESA_DIR=/Users/owner/Documents/Software/dev/test/mesa`,
`MESASDK_ROOT=/Applications/mesasdk` and `GYRE_DIR` unset.
The installed libstar.a members match the newly built hydro_rsp2.o and
set_flags.o byte for byte. See `install.log` and `install_check.json` in
the same review directory. The case was neither relinked nor rerun.
The user will relink; improvement over the saved 2000-step baseline
remains to be tested.

Both decay controls are 1d0 in the defaults and this case, as requested.
The RSP3 and star_LNA manuscripts now describe the conditional startup
seed consistently. Both PDFs were rebuilt; changed pages and adjacent
page breaks were rendered and inspected. RSP3 remains 12 pages and
star_LNA remains 35 pages.

## Crash reproduction after startup correction, 2026-09-19

The user requested another run to investigate a crash near models 1900--2000.
The current case executable was relinked by the user after installation.
Run a copy in `output/review/rsp3_crash_20260919/run` with eight threads and
the explicit checkout MESA_DIR. Keep physics and solver controls unchanged;
set only max_model_number=2300, profile_interval=1 and the Pi/Phi/covariance
profile columns. Preserve all original case inputs and outputs. No compilation
or relinking is needed. The manifest records the executable hash and inputs.

- [x] Reproduce the reported failure or reach model 2300.
- [x] Align retry reasons, residual rows, timesteps and face moment profiles.
- [x] Trace the failing terms and distinguish an implementation error from
      a limitation of the chosen closure or nonlinear solve.

The full result is `output/review/rsp3_crash_20260919/analysis.md`.
Model 1975 is the last accepted state; retries on model 1976 reach
dt=8.638468273433028e-7 s and stop at min_timestep_limit. A restart from the
same run's model-1000 photo reproduces the accepted history exactly and the
same final failure. Original case files and outputs remain intact.

The failing cell 75 has w=1.528268638841817e-23 cm/s while both bounding face
energies remain finite. Pi and Phi are finite and Phi is positive. The
assembled energy residual multiplied by csound/w loses its storage
derivative through cancellation of the explicitly linear source terms.
The native partial check on the saved model-1975 photo reports a diagonal
-7.739600595055e-5. A standalone reconstruction matches it, while the same
equation with common w factors canceled and 90-digit arithmetic gives
+2.3944755541716443e-7. The wrong sign is a confirmed numerical defect.
The numerical finite-difference column is unresolved at this tiny w and
is not used as the reference.

The next correction should cancel the known w factors analytically before
differentiation, reusing Source_div_w, D_div_w and v-grid Eq_div_w. Preserve
the original energy row, accepted state and time weights. Nonlocal transport
and u-grid heating cannot be divided analytically by w without a separate
derivation. No source correction, new floor, compile or relink was made
during this run audit. A corrected run is needed to establish that the
whole collapse is resolved.

## Small-w residual correction, 2026-09-19

The user authorized the correction and manuscript update. Preserve the
original energy equation, time weights, cell and face locations, zero-state
selection and positive roots. For 0<w<csound, construct R_energy/w term by
term before AD differentiation. Cancel the explicit w in Source, D, Dr,
current turbulent pressure and v-grid Eq. Retain accepted energy and pressure
history divided by current w. For u-grid Eq and current Lt, expand each
face velocity ratio as its local weight plus the neighbor contribution
divided by cell w. This retains transport and neighbor derivatives without
numerically canceling the local velocity factor. Accepted Lt remains an
independent history term. Keep the existing raw row at w>=csound and the
existing zero-state rules.

- [x] Implement the factorization without changing the physical fluxes or Eq/Uq.
- [x] Check values and derivatives with time centering, both velocity grids,
      nonzero transport, source seed and both RSP2 closures.
- [x] Verify the diagnosed tiny-w diagonal with the equivalent factored row.
- [x] Update and visually check both manuscripts, showing the cancellation.

Before-state files are in `output/review/rsp3_small_w_fix_20260919/before`.

The change is confined to `star/private/hydro_rsp2.f90`:
`do1_turbulent_energy_eqn` assembles the divided terms; its existing
`setup_dt_dLt_dm_ad` retains the luminosity time weight. Optional `k_div_w`
in `compute_Lt`, `compute_Eq_face` and `compute_Chi_face` expands the face
velocity ratio before differentiating it. Normal calls still return physical
luminosity, stress and heating. Divided calls do not store divided luminosity
in the physical Lt cache. `compute_Eq_div_w_cell` covers both velocity grids.

The explicit cancellations are

```math
(w^2-w_start^2)/w = w-w_start^2/w,
d/dw [w-w_start^2/w] = 1+(w_start/w)^2,

Source_3eq,k/w_k = (1/2) sum_f [-(d_r P)_f/rho_f]
                  [chi_T/(chi_rho Cp)]_f Pi_f/sqrt(e_t,f),
D_k/w_k = C_D w_k^2/Lambda_k,
Dr_k/w_k = w_k/tau_rad,k,

[a_f*w_k+(1-a_f)*w_(k-1)]/w_k = a_f+(1-a_f)*w_(k-1)/w_k.
```

Dr is omitted in the three equation mode. Only the explicit source numerator
w cancels; the face energy and its AD derivatives remain. Current pressure
divided by w is `(2/3)*alfap*rho*w`. Accepted pressure and energy still divide
by current w. Face ratios retain the neighbor contribution in current Lt and
u-grid stress. At the excised inner stress boundary, the ratio w_nz/w_nz is
exactly one. No floor, new dormant branch, physical closure, time weight or
continuous LNA equation is introduced.

`notes/check_rsp2_factored_energy.py` checks 864 residual comparisons and nine
partials per case with both closures, both velocity grids, zero/nonzero alfat
and alfam, source seed, pressure and luminosity time weights, and velocity
centering. Maximum scaled derivative difference is 6.11e-16. The failure
fixture gives +2.3944755541716443e-7, matching the independently derived
90-digit reference; the prior native derivative was -7.739600595055e-5.
Results are in `factored_energy_checks.json` in the review directory.
The existing three equation standalone suite and Fortitude pass. These checks
do not execute the corrected Fortran AD or establish full MESA convergence.

Installation succeeded with the checkout MESA_DIR, MESASDK_ROOT=/Applications/mesasdk,
GYRE_DIR unset and NPROCS=10. `install.log` and `install_check.json` record it.
The rebuilt hydro_rsp2.o matches its member in `build/star/lib/libstar.a`
byte for byte. The user's case was not relinked or rerun; the standing
instruction leaves relinking to the user. A corrected run past model 1976
remains necessary before claiming that the full timestep collapse is fixed.

The RSP3 PDF is 13 pages. Section 9, pages 9--10, shows the storage, source,
sink, pressure, face velocity and luminosity history cancellations and the
complete divided row. The star_LNA PDF remains 35 pages; page 31 states the
same evaluation and distinguishes it from the continuous LNA operator.
Both were rebuilt with settled references. Changed pages and adjacent pages
were rendered and inspected; no overflow or stray continuation page remains.

## Small-w AD overflow after factorization, 2026-09-19

The user's pasted run relinks the corrected library, then fails at model 47
with NaN derivatives in the cell-75 turbulent energy row. The accepted
history ends at model 46. This supersedes the prior pending-runtime status;
the factorization alone is not a complete correction.

The native star AD quotient forms `1/pow2(denominator)` in its derivatives.
This can overflow even when the quotient and its mathematical derivatives
are finite. The new row still uses that operation for accepted energy
divided by current w. Investigate with an isolated copy of the user's
executable, retaining the current ten-thread configuration, physics and
solver controls. Save profiles/photos every model; change only the stop
to model 100 and diagnostic output. Preserve the original case.

- [x] Reproduce the early failure and measure the cell-75 w trajectory.
- [x] Remove the avoidable intermediate overflow without a velocity floor,
      retaining accepted energy and all required derivatives.
- [x] Check native AD at extreme w, then install with the checkout MESA_DIR.
- [x] Update the manuscript and distinguish checked algebra from runtime evidence.

Evidence and before-state source are in
`output/review/rsp3_small_w_nan_20260919`.

The isolated run matches every accepted history column of the user's run
exactly through model 46, including the same final failure. Cell 75 declines
from w=6.78369e-10 at model 1 to 6.78369e-152 at model 45 and 6.78369e-154
at model 46. Its next downward trial overflows the native quotient derivative.
The earlier standalone checks did not extend to this range and missed it.

Keep the analytic factorization and evaluate remaining quotients using

```math
q = x/w,
dq = (dx-q*dw)/w.
```

`hydro_rsp2:div_by_w` implements this identity directly on the AD value and
partial array. It never forms an inverse squared denominator. It covers
accepted energy/pressure, accepted Lt, the source seed, both face velocity
ratios, and the Pi/w_face and w_cell/w_face source reconstructions. No generic
AD library or physical coefficient is changed. True unrepresentable quotients
or derivatives are not made finite by this identity.

Use `get_etrb_start(s,k)==0` for the local dormant branch and the raw zero-state
test. It is the accepted energy used in the original residual. Positive
representable energy still requires the positive branch; when w_start squared
underflows to exactly zero, the existing dormant row can reach w=0. This is
machine underflow of energy, not a configurable velocity floor. Nonzero Lt
and u-grid heating retain the existing restrictions on local branch selection.

`native_ad_check.f90` contains the helper extracted from the actual source and
uses the compiled native `auto_diff_real_star_order1` type and operators.
It reproduces five nonfinite legacy cases and zero failures for the new
division. Tests include zero numerators, equal tiny velocities, independent
partials, the model-46 history term, and the energy-underflow branch. At the
next model-46 trial its finite history partial is -100. The model-1975 storage
diagonal remains +2.394475554171644e-7. Tested velocities extend to 1e-310 for
representable quotient values and derivatives. This standalone native test
does not execute the complete corrected stellar model.

Installation succeeded with MESA_DIR=/Users/owner/Documents/Software/dev/test/mesa,
MESASDK_ROOT=/Applications/mesasdk and GYRE_DIR unset. The installed archive's
hydro_rsp2.o matches the rebuilt object exactly. Fortitude, the existing
three equation standalone suite, the 864 factored-residual checks and
`git diff --check` pass. Original case input hashes are unchanged.
The corrected stellar model was not run; the user retains case relinking.

The user's trace contains 14 bad-entropy-moment retries before the persistent
NaN failure and 22 energy-row-error retries. The latter failure is reproduced
and its arithmetic defect is addressed. The former retries are a separate
remaining diagnostic question, not demonstrated fixed by the quotient test.

Both manuscripts now document the quotient derivative ordering and exact
energy-underflow boundary. The rebuilt PDFs remain 13 and 35 pages. RSP3
pages 9--11 and star_LNA pages 31--32 were rendered and visually checked;
neither build has overflow or unresolved references.

## LNA consistency audit after the small-w corrections, 2026-09-19

The continuous positive-energy equations remain shared and consistent:

- `star_LNA_turbulence_closures:rsp2_source_for_star_LNA` calls the same
  `hydro_rsp2:compute_Source`; its moment reconstruction now uses `div_by_w`.
- D and Dr also call their nonlinear physical helpers. Lt uses the physical
  cached AD flux; optional divided evaluations do not replace that cache.
- `star_LNA_support:assemble_rsp2_moment_rows` calls `rsp2_moment_rhs` with
  instantaneous gradients and unit Pi/Phi inertia. There is no finite-dt
  history or solver residual normalization in the continuous operator.
- At a stationary positive background, canceling/dividing the energy row is
  a nonsingular row scaling of both stiffness and inertia, preserving its
  eigenvalues. This statement requires an equilibrium background.

The zero-energy treatment is not yet completely aligned. Nonlinear RSP2
uses zero represented accepted energy, including underflow, to permit its
local dormant branch. `rsp2_zero_w_for_star_LNA` still requires exact w=0
and both bounding moment faces dormant, apart from forced cells. Its
unconstrained row retains `delta(w*w)=2*w*delta w`.

A concrete local counterexample has w=1e-170, no transport, no source seed,
face-velocity hydro, finite neighboring face energies and Source/w=-1e-3.
The stored cell energy is exactly zero in double precision. Nonlinear
RSP2 selects its w=0 branch; LNA retains inertia 2e-170 and a spurious
local decay eigenvalue about -5e166 s^-1. This is not just row scaling.
At exact w=0 with nonzero Source/w, the raw LNA row can still be equivalent
to delta w=0; that does not cover a vanishing source coefficient or the
energy-underflow example.

Before claiming full boundary consistency, align LNA's dormant selection
with the local nonlinear branch while retaining active neighboring moment
rows and nonlocal transport. Do not copy a timestep predictor or blindly
freeze every zero-w cell into the continuous eigenproblem. This audit makes
no further source change and does not claim a new eigenmode validation.

## LNA and remesh boundary corrections, 2026-09-19

Implementation plan before editing:

1. Use current represented cell energy, `get_etrb(s,k) = w(k)**2`, in
   `star_LNA_support:rsp2_zero_w_for_star_LNA`. Preserve the fully dormant
   three equation branch. For a local row without a source seed or kinetic
   transport, select `delta w = 0` when `Source/w <= 0` at zero energy.
   At a static background, viscous heating and pressure work do not supply
   a finite source there. Positive driving and imported energy must retain
   the physical row. Propagate errors through the equation map, assembly
   and equilibrium diagnostic. Do not use accepted timestep history in LNA.
2. Make the shared dormant moment test use represented adjacent energies,
   including underflow. Pi and Phi must still both be zero. An active
   moment face next to finite energy remains dynamic.
3. Consolidate the post-remesh moment boundary checks in
   `hydro_rsp2:remesh_rsp2_moments`, using explicit new `nz`, `dq` and `xh`.
   All three callers must remap cell energy before this check. With the
   same face weights as the physical source,

   ```math
   e_{t,f}=a_f w_k^2+(1-a_f)w_{k-1}^2,\qquad
   e_{t,f}=0\ \Longrightarrow\ \mathrm{Pi}_f=0.
   ```

   This is the zero kinetic energy endpoint, not a finite energy flux
   limiter. Keep interpolated Phi on such an unforced face, so buoyancy
   can restart convection. Reapply forced masks and reject negative Phi.
   Do not change the conservative cell energy remap or interpolate face
   moments as cell integrals.
4. Verify a split of a quiet cell bordering turbulent cells, signed Pi,
   surviving Phi, unequal mass weights, underflow, and all three call
   paths. Check the LNA constraint and active branch separately, including
   transport and current versus accepted energy. Run existing regressions,
   compile/install with the correct MESA_DIR, and update both manuscripts.
   Do not relink the user's case or claim a full model or eigenmode test.

Completed implementation:

- `star_LNA_support:rsp2_zero_w_for_star_LNA` now uses current represented
  energy and the shared factored source. The equation map, matrix assembly
  and equilibrium diagnostic use that same decision and propagate `ierr`.
  Both RSP2 closures and velocity grids are covered. Mean-flow viscous
  heating is quadratic on the static LNA background, so it supplies no
  finite driving at this endpoint. The fully dormant three equation
  branch remains available with kinetic transport enabled when its
  bounding faces are also dormant.
- `assemble_rsp2_moment_rows` keeps active faces at positive face energy.
  It explicitly rejects active moments at zero face energy, where the
  square-root moment coefficients and reconstructed source do not have
  a regular linearization. This avoids silently freezing an active Phi
  reservoir into a spurious static mode. Such a state requires nonlinear
  evolution before a regular LNA background can be formed.
- `hydro_rsp2:rsp2_dormant_moments` tests represented adjacent energies.
  Nonzero Pi or Phi still prevents the dormant moment branch.
- `hydro_rsp2:remesh_rsp2_moments` owns the shared post-remesh check.
  `mesh_adjust:do_mesh_adjust` passes the new mesh arrays after remapping
  w, Pi and Phi. `tdc_hydro_support:remesh_for_TDC_pulsations` calls it
  after conservative w remapping. `adjust_mesh_split_merge:remesh_split_merge`
  calls it after AMR and synchronizes Pi/Phi from xh before derived state
  is rebuilt. The duplicated forced-face checks were removed.
- No new control, finite-energy limiter, change to energy conservation,
  interpolation order, or timestep centering was introduced.

Validation and installation:

- `output/review/rsp3_lna_remesh_20260919/native_check.txt`: 87 checks pass
  with bounds checking. Remesh and dormant-face tests call the installed
  MESA routines and native monotonic cubic interpolation. They reproduce
  the invalid flux inside a split zero-energy cell, retain Phi, preserve
  signed active fluxes and integrated cell energy, exercise unequal
  weights, merging, underflow, forced boundaries and negative variance.
  Stale old `s%nz` and physical arrays are deliberately supplied to verify
  that the remesh helper uses its explicit new arrays.
- The private LNA selector is extracted verbatim for the native harness;
  its Source/w evaluator is replaced by prescribed coefficients. This
  isolates classification, source-error propagation, independence from
  accepted history, both closures/grids and the transport exception.
  The earlier underflow counterexample changes from a spurious finite
  decay eigenvalue near -5e166 to an algebraic zero-inertia row. This is
  a branch test, not a coupled MESA eigenmode calculation.
- `notes/check_rsp2_remesh.py`, `notes/check_rsp2_three_equation_implementation.py`
  and `notes/check_rsp2_factored_energy.py` pass. The latter retains all
  864 residual/derivative comparisons; module dependency checks find no
  cycle. Fortitude passes for the five changed Fortran files, and
  `git diff --check` is clean. `source_changes.diff` in the review directory
  isolates this correction from the earlier uncommitted implementation.
- Full `./install` passed with
  `MESA_DIR=/Users/owner/Documents/Software/dev/test/mesa`,
  `MESASDK_ROOT=/Applications/mesasdk`, GYRE_DIR unset and NPROCS=10.
  The installed library is `build/star/lib/libstar.a`. Installation's
  standard package checks passed. No user case was relinked or evolved.
- Both manuscripts now state the dormant LNA selection and its physical
  limits. RSP3 also gives the common mesh endpoint rule. PDFs remain
  13 and 35 pages. RSP3 pages 10--13 and star_LNA pages 13--15 and 30--32
  were rendered and visually checked; neither build reports overflow or
  unresolved references. Full corrected stellar evolution and eigenmode
  validation remain pending.

## Ordinary remesh startup failure, 2026-09-19

The user's 22:34 run accepts model 1, then stops with `nonzero_ierr` before
the second solve. An isolated copy of the user's freshly linked executable,
with `report_ierr` enabled, reproduces all five first-step iterations and
shows `unpack_xh` failing after ordinary remesh. No user input or executable
was changed. The envelope remesh is scheduled for model 100, so it has not
run at this failure.

Reading the diagnostic model-1 photo gives Phi=1.9482144454916305e-11 at
old face 76 and Phi=1.5841980239542095e7 at face 75. Native MESA monotonic
cubic interpolation at the exact coordinate of face 76 returns
-1.862645149230957e-9. The polynomial evaluation loses the small positive
endpoint by cancellation; the accepted moment itself is valid. The earlier
small remesh fixtures did not cover this dynamic range.

`mesh_adjust:failed` calls `dealloc`, which passes the caller's `ierr` to
`do_work_arrays`. The successful deallocation resets it to zero and hides
the negative-Phi error. `unpack_xh` later rejects the resulting state.

Plan before editing:

1. Share the normal and envelope RSP2 face interpolation in
   `hydro_rsp2:interpolate_rsp2_face`, using native `interpolate_vector_pm`.
   Copy an old face value exactly when the requested coordinate is that
   face. For Phi between old faces enforce the monotonic interpolant's
   existing endpoint bounds:

   ```math
   \min(\mathrm{Phi}_j,\mathrm{Phi}_{j+1})\leq
   \mathrm{Phi}(q)\leq\max(\mathrm{Phi}_j,\mathrm{Phi}_{j+1}).
   ```

   Reject invalid old Phi; do not impose an arbitrary positive floor.
   Keep signed Pi and Y and the conservative w-squared remap unchanged.
2. Preserve the original error through `mesh_adjust:dealloc`, and let
   `report_ierr` identify the failing remesh operation.
3. Add a native regression from the exact model-1 values and all original
   face coordinates, plus points between faces and invalid input. Repeat
   existing checks and install with the correct MESA_DIR. Do not relink
   the user's case; distinguish native interpolation validation from a
   complete corrected run.

Completed:

- `interpolate_rsp2_face` now wraps the existing MESA monotonic cubic and
  restores exact old face values for Y, Pi and Phi. Between faces only Phi
  is bounded by its donor endpoint values. Invalid old Phi, including NaN,
  is rejected. The ordinary and envelope remesh call the same routine.
  AMR's convex Phi interpolation is unchanged.
- `mesh_adjust:dealloc` preserves an existing error while still returning
  a cleanup error when there was no earlier failure. The `failed` helper
  now honors `report_ierr` when identifying the remesh operation.
- The native regression reads the actual model-1 photo data and samples
  14,901 coordinates. The legacy interpolation reproduces negative Phi;
  the corrected value at face 76 is exactly the stored positive value.
  All 150 old face values, including the added inner endpoint, are exact
  for Y, Pi and Phi. Phi obeys the local endpoint bounds; between old
  faces, signed Pi and Y are identical to native cubic interpolation.
  Tests also cover dimensional envelope mass coordinates, all-zero Phi,
  negative/NaN input, and error preservation during cleanup. Results are
  in `output/review/rsp3_remesh_startup_20260919/native_check.txt`.
- Existing remesh and three equation regression scripts pass. Fortitude
  and `git diff --check` pass. Full installation passed with the explicit
  test MESA_DIR, MESASDK_ROOT, GYRE_DIR unset and NPROCS=10. The native
  interpolation test links against the newly installed MESA library.
- The original case inputs and executable match the recorded SHA256
  hashes. The isolated uncorrected executable reproduced the startup
  failure, but no corrected full case was relinked or evolved. The user
  will relink. This corrects numerical remap evaluation and error handling;
  it does not change the equations, time weights, turbulent energy remap
  or finite-energy convective flux.

## Authorized remesh evolution test, 2026-09-19

The user explicitly requested running, testing and fixing the case after
the interpolation correction. The case executable was relinked by the
user at 22:46, after the 22:43 installation. Start with that executable
in `output/review/rsp3_remesh_run_20260919/before`; keep the user's inputs,
outputs and executable untouched. The copy uses report_ierr, frequent
photos/profiles, no plotting and a 2000-model test limit. Physical controls,
ordinary remeshing through model 100, envelope remeshing at 100, and LNA
at model 200 remain as supplied. Use eight threads for the evolution run.

Before any further source edit, identify the first remaining failure and
its state transition or residual from this run. Preserve a restart photo,
derive the correction against the existing equations and MESA data paths,
then install and validate the actual evolution through the failed point.
The current request authorizes building and running an isolated test case
against subsequent fixes. It does not authorize overwriting the user's
case or changing its physical controls to make the test pass.

The reproduced failure is model 2, immediately after 149 cells become 350.
The first Newton correction clips Phi to zero, but the scaled correction
round trip gives, for example, `-4.1359030627651384d-25` from an accepted
`3.4660054625927101d-9`. Later retries produce Pi of order `1d-11` on faces
with exactly zero accepted and current kinetic energy and moments. These
faces have algebraic zero moment rows. Neither error is repaired by
reducing the timestep. The user's independent run without remeshing
advances, so the remapped initial state also needs checking.

Next correction and checks:

1. In `solver_support:Bdomain`, evaluate the same update as
   `star_solver:apply_coeff` and `set_vars_for_solver`:
   `xh_start + (solver_dx + correction_factor*x_scale*B)`. Round a clipped
   correction toward the admissible side if this reconstruction still
   crosses the boundary. Do not add a physical Phi or w floor.
2. Preserve the exact zero correction required by a dormant zero moment
   row when both accepted moments are zero. This only removes linear solve
   roundoff; an active face retains its dynamic equations. New kinetic
   energy can activate the face on the next Newton evaluation.
3. Capture the remapped state before the predictor, check turbulent energy
   conservation and face interpolation, then rerun the same remeshing case.
   These fixes are not yet validated by evolution.

The state dump identifies a larger remap inconsistency at new face 113.
Both new adjacent cells have `w=6.78368885d-10`, but directly interpolated
Pi is `-1.70077113d8`. The old boundary face had a turbulent outer cell
(`w=2.02648373d5`). Inserting a face inside its nearly quiet inner cell
therefore imports a boundary flux without its kinetic energy. Monotonic
interpolation of each variable separately cannot prevent this.

Test a coupled remap using the existing mass or arithmetic face energy:

```math
w_f=\sqrt{a_f w_k^2+(1-a_f)w_{k-1}^2},\qquad
\mathrm{Pi}_{\rm new}=w_{f,\rm new}\,
 \mathcal I\left[\frac{\mathrm{Pi}_{\rm old}}{w_{f,\rm old}}\right].
```

Set the old ratio to zero where the represented old face energy is zero.
Keep the conservative cell energy remap and the independent Phi reservoir.
This transfers the correlated entropy amplitude and reconstructs Pi with
the kinetic energy actually present on the new mesh. It changes the mesh
transfer only, not the finite energy evolution equations or a flux limiter.
Validate the ordinary mesh first before extending to the other mesh paths.

The first coupled-remap trial removes the large Pi/w mismatch, but the
zero-energy trial-state failure remains. Testing the *previous* dormant
state in Bdomain was insufficient because the shear predictor seeds w
before the Newton solve. Remove that unsuccessful guard. After applying
a Newton correction, if both represented cell energies and both accepted
moments are exactly zero, impose the homogeneous local solution Pi=Phi=0
exactly in solver_dx. A face with accepted Phi retains its reservoir and
requires positive trial w. This follows the existing zero-moment equation
branch; it does not cap finite-energy Pi. Test this combination next.

The combined trial passes model 2 but stalls at model 22. The polynomial
bound correction needs to advance the *total solver displacement* to the
next representable value, not repeatedly advance B: after cancellation
against xh_start, advancing B may require an impractical number of ulps.
The revised Bdomain passes 3520 native reconstruction checks, including
accumulated corrections and smaller line-search coefficients.

The later failure is preceded by internal luminosities of about 1d6 Lsun
in the nominal 60 Lsun model. Model 10 has artificial convection in the
formerly quiet deep envelope. `mesh_adjust:do1_lnT` with
`mesh_adjust_get_T_from_E=.false.` calls `get_old_value_integral`, which
uses a constant value in each old cell. Every child initially gets its
parent's identical temperature. This creates temperature plateaus and
sharp gradients at the old cell edges. Banded and bcyclic restart tests
both fail at model 22; disabling dynamic gradL reaches 50 but still has
31 retries and a 2.29d-5 s timestep. Neither is a solution.

Next test: in the RSP2 non-energy temperature remap, use existing
`get1_lpp` and `get_xq_integral` to reconstruct and integrate lnT. Preserve
its old cell averages. At an excised inner boundary use the one-sided
last-cell slope instead of imposing a flat central cell. Keep the actual
center boundary, conservative energy option and cell w-squared remap
unchanged. This follows the documented polynomial temperature remap;
it does not impose hydrostatic equilibrium on an evolving pulsation.

The polynomial lnT run passes ordinary remeshing through model 100,
envelope remeshing at 100 and reaches LNA at 200. At model 10 its maximum
internal luminosity is 96.73 Lsun, compared with 1.04d6 Lsun with the
constant-cell temperature remap. The test continues to 2000.

The user asked whether the false temperature option should work. Yes;
it is a separate supported remap choice. A paired true-option test also
shows timestep loss because its energy remap uses constant old cell
values. Consolidate reconstruction before `do1_lnT`: integrate the same
native polynomial representation of lnT for false, or internal energy for
true. The latter preserves the mass integral of internal energy and then
uses the existing EOS inversion and optional PE/KE correction. Do not
change the chosen energy equation or any runtime physics controls.

The first polynomial temperature run completed 2000 models, with ordinary
remeshing through 100, envelope remeshing at 100, LNA at 200 and time
centering from 300. It ended at the 30.781 s timestep cap with 5 iterations
and 68 cumulative retries. This establishes recovery, not retry-free
convergence. Artifacts are in `output/review/rsp3_remesh_run_20260919/mesh_all`.

The direct reconstruction check found that `get1_lpp` can limit curvature
without recentering its constant coefficient. The RSP2 caller therefore
sets `c0 = old_average - c2*dq_old**2/24` after reconstruction. This is
required to retain the old cell integral when its limiter changes c2.
The native Fortran check integrates 40000 positive child averages on
unequal meshes; relative integral error is below 6.5d-14. Without
recentering, a sampled cell-average discrepancy reaches 4.15 percent.
The existing generic composition reconstruction is outside this change.
Repeat both full tests against this conservative reconstruction.

Both-mode thermal reconstruction with consistent cell averages:
`conservative_false` completed 2000 models at 30.7814 s, 68 retries and
20.3913 days. Its minimum timestep was 9.83368 s. `amr` completed 100
models with 2183 zones, 16 retries and a recovered 54.5907 s timestep.
The true-option run with constant old mechanical energy failed at 564;
it had only advanced to 0.142556 days. The isolated comparison
`energy_without_PE_KE`, changing only
`max_rel_delta_IE_for_mesh_total_energy_balance=0`, reached 100 at
34665.8 s and 20.6268 days with 13 retries. This isolates the mechanical
energy correction as the remaining true-option problem.

For RSP2, reconstruct and integrate both internal energy and old
specific PE+KE using the same cell-average polynomials. In `do1_lnT`, use

```math
u_{k,\mathrm{new}}=\overline{u}_{k,\mathrm{old}}+
 \operatorname{clip}\!\left[
 \overline{(e_{\mathrm{grav}}+e_{\mathrm{kin}})}_{k,\mathrm{old}}
 -(e_{\mathrm{grav},k,\mathrm{new}}+e_{\mathrm{kin},k,\mathrm{new}}),
 \pm f\overline{u}_{k,\mathrm{old}}
ight],
```

where f is the existing maximum relative energy correction control.
The overbar is the new-cell mass average of the reconstructed old state.
Before the existing cap and EOS fallback, this preserves total energy;
using constant old PE+KE creates artificial differential heating during
splits. Cell turbulent energy already has its own conservative remap
and is not counted again. Generic non-RSP2 behavior remains unchanged.
Test the corrected true option without disabling the energy correction.

The consistent mechanical reconstruction passes model 100 with the energy
correction enabled: dt=1d5 s, age=31.2828 days, 11 retries. It passes LNA
at 200 and the former model-564 failure, with a recovered timestep above
70 s and 43 retries near model 580. Continue to 2000.

Native checks against the installed RSP2 face routines reproduce a direct
Pi interpolation of -1.5d8 inside a cell with w=1d-9, compared with
-1.06066d-6 from the coupled transfer. Both face weighting choices pass,
cell turbulent energy is unchanged, and zero-energy Phi reservoirs are
retained. The linear thermal plus mechanical reconstruction has absolute
error 2.85d-14 in the analytic fixture. Native Bdomain checks cover 4400
scaled/accumulated corrections. These checks and the source regressions
are recorded under `output/review/rsp3_remesh_run_20260919`.

### Verified remesh result

Installed successfully with
`MESA_DIR=/Users/owner/Documents/Software/dev/test/mesa`; final build log
is `output/review/rsp3_remesh_run_20260919/install_mechanical_final.log`.
The user case was not relinked by this audit. All runs used isolated
copies and eight threads. The user's own inlist/executable changed during
the audit and were left alone.

| Thermal remap | Models | Minimum dt (s) | Final dt (s) | Retries | Median iterations |
| --- | ---: | ---: | ---: | ---: | ---: |
| lnT, control false | 2000 | 9.83368 | 30.7814 | 68 | 5 |
| Internal energy, control true | 2000 | 19.6903 | 28.3263 | 43 | 5 |

Both runs include ordinary remeshing through model 100, the one-time
350-zone envelope remesh, LNA at model 200 and the change to time
centering from model 300. The false run advances 20.3913 days; the true
run advances 33.0678 days. These are convergence checks at equal model
count, not a comparison at equal age or a validation of pulsation growth
rates. The new mechanical-energy branch is true-only and leaves the
verified false path unchanged. The 100-model split-and-merge test reaches
2183 zones, with both splits and merges, 16 retries and dt=54.5907 s.
That is a shorter AMR startup test, not a 2000-model AMR validation.

The complete numerical result is in `verified_results.json` beneath the
review directory. Source regressions, six-file Fortitude and
`git diff --check` pass. The RSP3 manuscript now describes the coupled
Pi remap and both thermal reconstructions; the rebuilt PDF is 14 pages
and its changed pages were checked visually. The continuous evolution
and LNA equations have not changed in this remesh correction.

Principal source locations:

- `star/private/mesh_adjust.f90`: `do_mesh_adjust`, `do1_lnT`,
  `do_RSP2_face_var`.
- `star/private/hydro_rsp2.f90`: `rsp2_remesh_w_face`,
  `interpolate_rsp2_face`, `remesh_rsp2_moments`.
- `star/private/tdc_hydro_support.f90`: envelope Pi/w transfer after
  conservative cell energy remapping.
- `star/private/adjust_mesh_split_merge.f90`: `do_split`, `do_merge`.
- `star/private/solver_support.f90`: `Bdomain`, representable lower bounds.
- `star/private/star_solver.f90`: exact homogeneous moment trial at zero
  represented energy and zero accepted moments.

Diagnostic photos were thinned after completion: startup, mesh/LNA/time-centering
checkpoints, failure boundaries and final photos were retained. Terminal
logs, profiles, history, source snapshots and test inputs remain.

### Pi, Phi and RSP2 mixing labels

The public face arrays are `s% Pi` and `s% Phi`. The profile columns `Pi`,
`Phi` and `Pi_covariance_excess` are registered in `star_profile_def` and
read in `profile_getval`. The user's current case list does not enable
these columns. Its existing `Lc` column calls `get_Lconv`, which returns
`s% Lc` for either RSP2 model. `compute_Lc_terms` sets

```math
L_c=4\pi r^2\rho_f T_f\Pi_f,\qquad F_c=\rho_f T_f\Pi_f.
```

The current label in `mix_info:set_mixing_info` uses `abs(Lt/L)` before
testing `abs(Lc/L)`. This can label an unstable convective face as
overshooting and misses stable turbulent faces when `alfat=0`. The sign
of Lt gives the transport direction, not the sign of its divergence or
the local stability.

Replace this classification for both RSP2 forms with:

- No mixing label if `D_mix <= 0` or the face convective velocity does not
  exceed `RSP2_min_conv_vel_for_mixing_type`, default 1 cm/s.
- Overshoot label for active mixing with `Y_face < 0`.
- Convective label for active mixing with `Y_face >= 0`.

Here `Y_face=gradT-gradL`, using the neutral gradient of the active
temperature equation. It reduces to `gradT-grada` when composition and
dynamical corrections are absent and the gradient definitions coincide.
The label describes turbulence in a locally stable layer; in a pulsating
model this can also be turbulence retained from an earlier unstable phase.
It does not establish that nonlocal transport is its source. No sign
condition on Pi, Lc or Lt is imposed. The threshold suppresses labels for
negligible turbulent seeds without setting their velocity or D_mix to zero.

Remove the two old luminosity threshold controls, as requested. Preserve the native surface
plotting convention that copies the second face's label to the first.
Do not change the physical fluxes, diffusion coefficient or TDC labels.
`mixing_type` also sets convection boundaries, so other explicitly enabled
boundary mixing prescriptions can respond to the revised classification.

Validation: the scoped source diff and all five control wiring locations
were checked. No references to the removed controls remain in tracked
source or the case inlists. Fortitude and `git diff --check` pass.
Installation completed successfully with
`MESA_DIR=/Users/owner/Documents/Software/dev/test/mesa`; the build log is
`output/review/rsp2_mixing_labels_20260919/install.log`. The user's case
was not relinked or run, and its profile list was left unchanged. No new
evolutionary comparison was made for this classification change.

### Restart at model 115000, 2026-09-20

The user authorized `./re x00115000` for about 1500 steps to diagnose a
timestep collapse and stalled Newton iterations, with profiles every step.
Use an isolated copy of the current case, executable and photo under
`output/review/rsp3_restart_115000_20260920`. The baseline preserves the
physics controls: RSP3, alfat=0, v flag, dynamical gradL, dPrad/dm,
time centered pressure and luminosity, no remeshing. Its stopping model is
116500. Record each accepted model, save profiles every step with Pi/Phi,
and photos every 50 steps. Disable only graphical output in the copy.
The copied `re` runs the frozen executable through a run-only Makefile.
The original case and executable remain unchanged.

Plan: reproduce the failure, align retries with the responsible residual
and correction rows and their zones, then audit their derivatives, domain
limits and timestep terms. Compare alfat=0.1 from the identical photo.
Do not infer a zoning cause or tune tolerances without this evidence.
Preserve the continuous equations and conservative shared fluxes in any
numerical correction. Validation and results will be appended here.

Both runs reached model 116500 on the unchanged 149-zone mesh. Baseline:
80 retries, minimum dt=12.6542 s, median four accepted iterations and mean
6.3033. Alfat=0.1: 87 retries, minimum dt=11.3893 s, median four and mean
6.1747. All retry residuals are dPhi_dt at face 74 or 75; accepted stalled
steps also involve face 76. The variance residual reconstructed from the
profile at model 115550 is 3.81799136e-5, matching the solver. Phi is clipped
to zero while positive Pi in a stable layer continues to demand a negative
variance derivative. Tolerance relaxation accepts some such steps at
iteration ten; attempts at larger dt retry. This is not evidence for
changing zoning, luminosity normalization or the negative-w treatment.

The user authorized deriving a physical closure correction. The linked
derivation retains the default local stationary solution while making
local buoyancy, decay and mean compression covariance consistent. A
homogeneous ODE confirms the old closure crosses Phi=0 from an admissible
state without any mesh or Newton solver. The corrected local ODE does not.
The existing staggered source and energy-only nonlocal transport still
require compatible changes before this can be a complete MESA correction.
No physics/source patch or installation was made for this investigation.
All test edits are confined to the isolated copies. A final hash comparison
found that the user case inlist and executable changed during the audit;
those external changes were left intact. The profile list still matches
the original snapshot. The manifest identifies the frozen executable and
inputs actually used for both comparisons.

### Moment placement audit, 2026-09-20

The user requested an audit of the conflicting face/cell recommendations
for pulsation, evolution and later nonlocal moments. The audit compares
four layouts and corrects both earlier categorical recommendations.
The simple cell proposal loses convective coupling to alternating mean
entropy when gradients and Pi are both averaged. Its local storage/decay
rows and adding Pi diffusion do not restore that response.

The preferred simple long-term design places e_t, Pi and Phi together on
faces, keeping the gas thermodynamics in cells. This combines direct
thermal coupling with local covariance and common moment transport.
Conservative storage and transport can be projected from the face control
volumes into cells on unequal masses. Pressure work, Eq/Uq, boundaries,
remeshing, restart conversion and LNA still require a complete derivation.
TDC's face velocity is a precedent, not an interchangeable hydro closure.

Standalone checks reproduce the extra thermal null mode of the simple
cell scheme, the compact face response, common diffusion/remap positivity,
and conservative unequal-mass face energy projection. A separate
Crank-Nicolson example demonstrates why common spatial transport alone
does not guarantee nonnegative variance at large timesteps. It does not
justify changing total luminosity time centering.

Details and source references are in `rsp3_layout_audit.md`; reproducible
checks and results are `check_layouts.py` and `layout_checks.json` in the
existing restart audit output directory. No production source, case inputs,
compilation or new MESA execution was involved. The layout recommendation
is not a demonstrated improvement of the full stellar solver and is not
implemented.

### Alfat restart and PGSTAR clarification, 2026-09-20

The user reported that alfat appeared active after restarting x00115000
with RSP2_alfat=0, based on the pale overshooting region in PGSTAR and the
remaining stiffness. The current inlist sets zero, the photo does not
restore RSP2_alfat, and the case extras do not override it. `compute_Lt`
returns zero and clears `s% Lt(k)` when alfat is zero; `set_etrb_start_vars`
recomputes Lt_start with the current control.

Saved profiles from the prior isolated zero-alfat comparison have Lt and
Lt_start exactly zero at models 115001 and 116500. The nonzero-alfat
comparison has nonzero Lt in those profiles. No new MESA run was needed.
The current case's saved profile is model 1, not this reported restart,
so it is not evidence about the user's current live calculation.

PGSTAR's Mixing panel uses `clr_overshoot=clr_Beige`. RSP2 assigns that
label when D_mix is positive, convective velocity exceeds the labeling
threshold, and Y_face is negative. Neither Lt nor alfat enters this label.
The photo retains w, Pi and Phi; setting alfat to zero removes spatial
turbulent-energy transport, not the existing turbulence or its local
production and decay. The previously diagnosed variance failure occurs
at alfat=0 as well. There is no demonstrated restart-control defect to
patch. Clarified the existing alfat documentation without changing physics.

The velocity panel plots conv_vel/csound. In the saved zero-alfat run at
model 116500, face 140 has conv_vel=2.56049 cm/s and conv_vel/csound=1.82223e-7,
while log_D_mix=10.2641. Thus a velocity that appears zero on this linear
plot can still give finite D_mix=conv_vel*mixing_length/3. The face velocity
uses both neighboring cell w values. If both are zero, RSP2's mixing
contribution is zero. This diagnostic check does not establish that the
current moment closure predicts a physically correct turbulence tail.

A subsequent live profile `LOGS/profile2.data`, model 120000, confirms
Lt=Lt_start=0 while 104 of 149 faces carry the overshoot classification.
The user then restarted again; the saved profile remains identifiable by
its model header and hash. Selected values are preserved in
`output/review/rsp3_restart_115000_20260920/mixing_budget_120000.json`.
At zone 100, SOURCE=936.736, Eq=7022.948 and DAMP=2211.797 erg/g/s. Net
turbulent energy production is positive. At zone 60 the stored turbulence
is decaying, since SOURCE+Eq-DAMP is negative. The stable-layer mixing
therefore includes viscous production with alfar=0 and alfam=0.25, as well
as finite-lived turbulent energy. The label is a stratification category,
not evidence that Lt imported the turbulence. A finite local Pi source
can also contribute. This confirms the supply terms without endorsing
the known deficient variance closure or changing its physics.

The user requested a temporary checkpoint of the current RSP3 state before
the proposed layout change. Include the implementation, current RSP3 test
inlist, related source notes and manuscripts, and small audit evidence.
Retain the known closure limitation explicitly. Raw stellar outputs and
unrelated GYRE research are not part of this source checkpoint.

### Turbulence tails and the changed mixing label, 2026-09-20

Checkpoint `6d6bc2848` preserves the current implementation before a new
moment layout. The user questioned whether recent remeshing left turbulent
energy in stable layers. The following checks distinguish the plotted
category, the existing turbulence, and the separate remesh failure.

The earlier `mix_info:set_mixing_info` rule required a nonzero Lt/L ratio
to assign `overshoot_mixing`. With alfat=0 it could never draw that category,
even when w and D_mix were finite. The new rule labels finite turbulence
with Y_face<0 as overshooting. Applying both rules to the same model 120000
profile changes the recomputed overshoot count from 0 to 103; the stored
profile has 104, including the surface label copied from its neighbor.
No physical array is changed by this comparison. This establishes a cause
of the newly visible white regions, not the validity of the underlying
closure. Mixing classifications also enter boundary and mesh decisions,
so the source change cannot be called purely cosmetic in every run.

Saved model 110000 already has substantial deep turbulence while those
faces are labeled `no_mixing`. Model 120000 has exactly the same cell
masses. At zone 120, w changes from 2124.535 to 2485.164 cm/s and
log10(D_mix) from 13.32951 to 13.39437. At zone 140 the corresponding
values are 2.133599 to 2.638337 cm/s and 10.20360 to 10.29577. These
quantities were present before the new category appeared.

Pi and Phi are absent from the user's profile columns, so their current
values cannot be recovered from those files. The earlier isolated
zero-alfat restart profiles do include them. At model 115001, face 120
has Pi=-3.354153e5 erg cm/(g K s) and Phi=3.176764e9 erg^2/(g^2 K^2).
At face 140 they are -0.4000762 and 1240.278 in the same respective units.
They are small relative to the convective peaks but are not zero. That
restart retains a fixed mesh through 116500. Lt is zero throughout both
saved endpoint profiles. The small audit evidence, file hashes and cell
mass comparisons are in `mixing_tail_comparison.json` beside the existing
restart results. The earlier labeling comparison is in
`mixing_label_comparison_120000.json`.

The ordinary `mesh_adjust:do_etrb` and envelope
`tdc_hydro_support:remap1_cell_average2` transfer cell turbulent energy by
mass overlap:

```math
 (w_j^{\rm new})^2 = \frac{1}{\Delta m_j^{\rm new}}
 \sum_i \Delta m_{ij}^{\rm overlap}(w_i^{\rm old})^2.
```

This can spread an existing nonzero cell average across an overlapping
new cell. With zero energy in every contributing old cell it returns zero.
Pi and Phi are transferred separately on faces, with Pi reconstructed
using the remapped face energy. Their plotted support need not exactly
coincide with the cell w support. This source audit does not exclude
diffusion from repeated remapping or certify moment realizability.

The fresh attached `./rn` output is a separate startup: model 100 has
478 cells, followed by the special envelope remesh. The printed surface
pressure changes from 65.46836 to 8642.595 dyn/cm^2 and opacity from
7.506293e-3 to 5.492591e-5 cm^2/g. The next model fails with a flux
residual of 45.73 and eventually reaches the minimum timestep. This
locates a failure immediately after envelope reconstruction, but does
not identify the turbulence tail's origin. No new MESA execution,
compilation or production source change was made for these checks.

The next screenshot shows model 115600 at age 0.5663381 yr. The saved
isolated baseline `LOGS/profile600.data` matches that model and age, with
Pi and Phi peaks of 1.036717e15 and 1.586204e18. At face 120
(logtau=6.310665), their magnitudes are only 1.158683e-9 and 2.035645e-9
of those peaks, and Lc/L=-1.718178e-5. Describing them as essentially zero
on the plotted linear axes is appropriate. Their nonzero absolute values
alone do not explain or validate the large mixing coefficient.

At the same location conv_vel=2168.838 cm/s and mixing_length=3.227758e10
cm, giving D_mix=conv_vel*mixing_length/3=2.33351e13 cm^2/s. The cell has
w=2327.582 cm/s. The current mixing prescription depends on this kinetic
energy through w, independently of Pi and Phi; negligible heat transport
therefore does not make the prescribed chemical diffusivity negligible.
This is a property of the implemented closure, not evidence that this
stable-layer chemical transport is physically justified. No state reset
on restart or remeshing origin is established by the screenshot.

### Stable layer energy budget and closure audit, 2026-09-20

The user requested a cause and a fix. `audit_stable_tail.py` integrates
all 1500 previously saved baseline profiles, with no new MESA execution.
It checks identical cell masses and zero Lt at every step, and verifies
the cell equation

```math
 w_k^2-w_{k,\rm start}^2
 =\Delta t\,(\mathrm{SOURCE}_k+E_{q,k}-\mathrm{DAMP}_k)
```

for the tested cells. Turbulent pressure and radiative damping are zero
for this run. The maximum per-step error is below 1e-10 of the local
stored kinetic energy in every checked case. Over 126172.0676 s the
time-averaged rates, in erg/g/s, are:

| Cell | SOURCE | Eq | DAMP |
| --- | ---: | ---: | ---: |
| 90 | -6091.7695 | 19154.3974 | 12318.8362 |
| 100 | -1041.9830 | 3358.3282 | 2110.9173 |
| 120 | -1.042212 | 3.938718 | 0.866810 |
| 140 | -1.196198e-6 | 3.984962e-6 | 1.450673e-9 |

All four cells remain stable throughout. Their sampled face Phi values
remain positive; this deep-tail budget is not the separate clipped-Phi
failure at faces 74--76. The integrated rates reproduce their net energy
growth. This establishes local shear production as the ongoing energy
supply over the audited interval. It does not establish when every part
of the initial tail first formed. The newly visible category is still
explained by the label change above.

At model 115600, reconstruction of the Phi equation in cells 90, 100,
120 and 140 gives raw residuals below 5e-17 of their stored Phi. The weak
moments there are resolved solutions, not evidence that normalization
has left them unconstrained. The source audit also finds no clipping of
negative Pi or negative buoyancy production in the RSP3 cell source.

The local stationary moment equations explain why three-equation braking
can be much weaker than the one-equation source. For stable entropy
stratification, fixed thermodynamics, no mean strain and no radiation,
eliminate Phi with `Phi_rhs=0`, then Pi with `Pi_rhs=0`. In the small-w
limit with a finite stable entropy gradient,

```math
 \mathrm{SOURCE}\simeq-\frac{\alpha_\Phi}{3\Lambda}w^3,
 \qquad \mathrm{DAMP}=\frac{C_D}{\Lambda}w^3,
 \qquad E_q\ \mathrel{\propto}\ w.
```

Here alpha_Phi is the source coefficient
`RSP2_alfa_phi*4*sqrt(2/3)` and C_D is
`RSP2_alfad*(8/3)*sqrt(2/3)`. This is a uniform local asymptotic result,
not a replacement for the time-dependent staggered stellar equations.
The negative one-equation buoyancy source is proportional to w instead.
Consequently retaining that model's dissipation closure when adding the
independent moments permits a finite shear-supported state in a stable
layer. Zero Lt alone does not forbid that solution.

Kupka, Ahlborn and Weiss (2022), sections 3.1 and 4.2, independently identify
weak buoyant braking and excessive mixing in the original three-equation
Kuhfuss closure. Their section 3.6 introduces enhanced stable-layer
dissipation. This is a relevant physical correction candidate, not proof
that their stellar-evolution calibration fixes this pulsating envelope.
The primary source is `notes/references/kupka_ahlborn_weiss_2022_three_equation_kuhfuss.txt`
and https://arxiv.org/abs/2207.12296.

Changing the mixing label or zeroing w wherever Pi looks small would hide
this solution without correcting its equations. The published stable-layer
dissipation change is a physical closure extension. The user has been
asked whether to add it as an optional control or retain the closure while
completing the separately planned moment-consistency correction. No such
physics change has been made on the strength of this audit alone.

### Kovacs moment and enthalpy audit, 2026-09-20

The user requested a comparison before deciding on the Ahlborn stable
layer dissipation change. The equations, source mapping and independent
thermodynamic check are recorded in
[rsp3_kovacs_equation_audit.md](rsp3_kovacs_equation_audit.md).
The current local moment decay represents a closure of transport
divergence. It is distinct from the extra viscous moment decay in that
paper, and it is not a spatial variance flux that can be substituted
into an enthalpy flux formula. Saying that the uncomputed spatial flux
must be zero was too strong. Printed normalization inconsistencies must
also be resolved before adopting any of those extensions.

No source, controls, manuscript equations or case inputs were changed
for this comparison. No MESA compilation or execution was performed.

The follow-up comparison records that local transport relaxation and
explicit viscous moment decay have identical algebraic dependence here.
Adding the quoted viscous coefficients changes the default turnover
decay by about one percent; substituting them would greatly reduce it.
No replacement or additional damping was implemented.

For the user's requested explanation to revisit later, see
[Local relaxation versus viscous decay](rsp3_kovacs_equation_audit.md#local-relaxation-versus-viscous-decay).
That section records the plain language distinction, combined equations,
coefficient comparison, source locations and unchanged implementation status.

### Buoyancy dissipation proposal, 2026-09-20

The user asked about the published beta correction without moving the
moments. See [rsp3_buoyancy_dissipation_plan.md](rsp3_buoyancy_dissipation_plan.md).
The plan derives a noniterative dissipation product with finite first
derivatives at zero w, using the existing harmonic length and one new
c4 coefficient. It distinguishes the proposed dissipation-only scope
from shortening the common length in every closure. It also records
the start-of-step Brunt cache, composition gate, mesh lifecycle and LNA
decisions that remain part of the implementation. No physics was changed.

### 350-cell envelope remesh crash, 2026-09-20

The current case fails at model 101 immediately after the special envelope
remesh from 478 to 350 cells. An isolated restart from the same model-100
photo reaches model 150 when only that remesh is disabled. An earlier
executable also fails with the current photo. The failure is therefore not
caused by the latest mixing-label change.

The old cell 10 contains `w = 1.576520312470981e-141 cm/s`. The finer outer
grid puts new cells 72 and 73 wholly inside it. Conservative overlap preserves
its microscopic kinetic energy, while interpolated Phi is nonzero. Pi is
already remapped as Pi/w followed by multiplication by the new face w, so its
initial ratio is finite. In `rsp2_moment_source`, however, the divided energy
row has a derivative with respect to Pi proportional to inverse face w.
The Pi equation can create finite Pi through buoyancy times Phi.

`RSP2_adjust_vars_before_call_solver` supplied a finite trial state for that
case only when both adjacent w were exactly zero. In diagnostic copies,
changing only old cell 10 w to zero or to 1e-12 lets the 350-cell run reach
model 150; the unchanged photo fails at 101 with dt below 1e-6 s. The same
eight-thread setting was used for this comparison. These are diagnostic
photo edits, not a proposed reset of the user's model.

The small source correction extends the existing trial guess to w below
roundoff relative to the buoyancy impulse. With the face acceleration per
unit entropy returned by `rsp2_buoyancy_face`,

```math
 w_{\rm trial}=\Delta t\,|\texttt{rsp2_buoyancy_face}|\sqrt{\Phi},
 \qquad \max(w_k,w_{k-1})\le\epsilon_{\rm machine}w_{\rm trial}.
```

The existing predictor then sets the two trial velocities and increments
trial Pi by `dt*rsp2_buoyancy_face*Phi`. The frozen start state and all
residual equations are unchanged. There is no physical turbulence floor,
new closure, altered time centering or LNA equation change.

Status: the source patch is prepared and statically checked. It has not been
compiled, installed or tested with the unchanged photo. Compilation approval
is pending. The run matrix and isolated diagnostic build are recorded in
`output/review/rsp3_remesh_350_20260920/README.md`.

### RSP2 behavior with the three equation flag off, 2026-09-20

The user's new run reaches model 200 with dt 0.0487994191 s, 63 retries and
only 0.0184095523 days of evolution. An isolated run using an existing
diagnostic executable reproduces every accepted model exactly in the
principal history columns. The run and source audit are recorded in
`output/review/rsp2_off_regression_20260920/README.md`.

The positive-w divided residual and positive trial restriction introduced
during RSP3 development also applied to ordinary RSP2. At cell 310 the
original normalized energy residual is 3.7902e-113, but multiplication by
csound/w = 3.5206e107 gives the printed stalled residual 1.33437e-5.
The prepared source correction confines these two changes to the three
equation flag. Ordinary RSP2 retains the prior nonzero-energy residual,
the shared zero-state branch selection and nonnegative corrections.
Physical equations, Eq/Uq and time weights are unchanged.

The saved full LNA matrix also identifies the mode failure. The 0.702356-day
root is present, but its tiny w components are poorly resolved. Its normwise
error is 1.9548e-6, just outside the old refinement guard. A fixed-frequency
inverse iteration using the existing banded factorization repairs the
vector before allowing a Newton frequency update:

```math
 (A-\sigma B)z=Bx,\qquad x_{\rm trial}=z/z_j.
```

Use a weighted average of the old and repaired vectors for backtracking so
tiny components are not lost by subtraction. Retain the original-matrix
residual check and final user tolerance. The source also enforces the
mathematical upper bound of one on the finite componentwise residual.
This prevents roundoff-only rejection at a limit of 1d0 without treating
that limit as evidence of an accurate eigenvector.

On the saved MESA matrices, an independent SciPy/LAPACK check recovers all
15 selected modes in one fixed-frequency iteration each. The first residual
falls to 7.27e-10 and the largest selected residual is 3.12e-7. No frequency
changes in those repairs. This verifies the proposed matrix operation;
the edited Fortran still requires compilation and native validation.
The existing diagnostic executable reproduces the history and turbulent
rows, but its full LNA matrix is not bitwise identical to the current case
executable. Some momentum and temperature-gradient coefficients differ
slightly. These checks therefore establish the repair on the saved operator,
not a native test of the edited source. The existing factored-energy suite
also passes 864 cases with nine derivatives each; its maximum scaled
derivative difference is 6.10e-16.
Both PDFs were updated and their changed pages visually checked. No MESA
compilation or installation was performed; approval remains pending.

### Installation, 2026-09-20

Following the user's explicit installation instruction, `./install` completed
successfully with MESA_DIR=/Users/owner/Documents/Software/dev/test/mesa and
MESASDK_ROOT=/Applications/mesasdk. The build includes the RSP2 residual and
correction guards, the RSP3 remesh trial guess and the LNA vector refinement.
The standard star, astero and binary checks passed. The case was not relinked
or evolved during this installation; the corrected RSP2 and RSP3 evolutionary
behavior still requires verification. The log is
`output/review/rsp2_off_regression_20260920/install_20260920.log`.

### RSP3 modes and radius profiles, 2026-09-20

The installed executable reproduces the user's model-200 history and term
audit in an isolated replay. Expanding the selected list to 60 identifies
the fundamental at terminal mode 18, output/kick index 19, period
0.7016683444 days and residual 1.18e-11. Its nodeless displacement closely
matches the earlier refined 149-zone fundamental. The first 15 selected
roots exclude it. Other roots need physical classification beyond their
matrix residuals. The attempted frozen comparison stops at the existing
RSP2 guard and is not evidence for the identification.

An isolated copy of photo 1000 advances one accepted step with no new
retries and writes model 1001 with native Pi/Phi profile columns. Lt is
zero throughout, while the overshoot-labelled layers carry at most
|Lc/L|=1.845e-4. D_mix is proportional to velocity times mixing length.
The case has mix_factor=0, so chemical mixing is disabled.

The tail is not solely an absence of buoyancy braking. At r=3.9784 Rsun,
Source=-5.02e-4, Eq=1.32e-3 and DAMP=1.10e-8 erg/g/s. Their sum agrees
with the accepted turbulent-energy increase. The implemented eddy-viscous
term supplies this part of the tail. Extra stable-layer dissipation remains
a proposal, not an established cure. Figures, measurements and limitations
are in `output/review/rsp3_modes_profile_20260920/README.md`. The user's
case and source were not modified for this analysis.

The follow-up alternating-flux audit finds a specific spatial weakness:
cell energy averages its two bounding face Pi contributions, while face
Pi/Phi respond independently. In a uniform background an alternating Pi
perturbation cancels from the cell energy source. The mechanism is present
in the saved deep stable layers: 24 sign changes over 37 faces, with
99.9799 percent cancellation of the two source contributions at cell 295.
The actual source reconstruction agrees to 6.45e-16. This supports testing
colocated face moments before attributing the pattern solely to missing
stable-layer dissipation. Full evidence and limitations are appended to
`notes/rsp3_layout_audit.md`; no production source was changed.

The radiative cooling follow-up checks `rsp2_moment_rhs`: alfar=1 gives
gamma_r=2 sqrt(3), with losses -Pi/tau_rad and -2 Phi/tau_rad. The current
case leaves alfar at zero. On the saved model 1001, alfar=1 would give a
Phi radiative sink time of about 2148 days at 1.25 Rsun, 7.81 days at
3.98 Rsun, and 0.021 seconds near the surface. These estimates motivate
a cooling-enabled physical baseline, but do not establish a cure for the
alternating flux or its cell/face coupling. Details are in the audit README;
no new MESA run, control change or source change was made.

### Face w conversion plan, 2026-09-20

`notes/rsp2_face_w_implementation.md` expands the layout audit into a full
source-mapped plan. It specifies face residuals, cell energy quadrature,
dual-volume transport and the compatible gas-face Lt. It records the
remaining local pressure/viscous work allocation and boundary decisions,
all three remesh paths, versioned restart conversion, zero-state handling,
temperature-row coverage and LNA changes.

The eventual target is common face energy machinery for RSP2 and RSP3,
with separate closure equations and separate validation against current
cell RSP2. At the user's request the viscosity section targets TDC's
spatial discretization while retaining separate RSP2 routines. It identifies
which cell/face w ratios can disappear and preserves the distinctions in
normalization, time states, boundary stresses and mass corrections.
This is a documented proposal, not a production-source change.

`notes/check_rsp2_face_w_plan.py` verifies the candidate fixed-mass storage,
Lt, Eq, pressure allocation and common moment-overlap identities on unequal
cell masses. All checks pass; results are saved under
`output/review/rsp3_face_w_plan_20260920`. The check also demonstrates that
global pressure-work conservation does not imply local equality. Boundary
masks, full hydro work, AD and stellar validation remain to be completed.
