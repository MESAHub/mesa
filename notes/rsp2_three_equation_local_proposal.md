# RSP2: two additional local moment equations

Date: 2026-09-19. Status: literature check and proposed equations only.
No MESA source, controls, solver state, or manuscript equations are changed by
this proposal. No MESA compilation, installation, relinking, or model run.

Implementation follow-up: the user subsequently requested implementation.
The source audit, provisional edits, unresolved discretization questions and
ordered checklist now live in
[rsp2_three_equation_implementation.md](rsp2_three_equation_implementation.md).
This proposal remains the derivation record, not an implementation-status claim.

## Sources and scope

- Flaskamp (2003), *Nichtlokale und zeitabhängige Konvektion in Sternen*,
  printed pp. 22--28, especially equations (3.12)--(3.18), (3.21), and section
  3.6 on initialization and branches:
  <https://mediatum.ub.tum.de/doc/602965/602965.pdf>.
- Braun, Ahlborn, Kupka & Weiss (2026), section 2 and section 4.3:
  <https://arxiv.org/html/2604.06151v1>.
- Kupka, Ahlborn & Weiss (2022), appendix A:
  <https://doi.org/10.1051/0004-6361/202243125>.
  Local copy: `notes/references/kupka_ahlborn_weiss_2022_three_equation_kuhfuss.pdf`.
- Ahlborn et al. (2026), equations (5)--(10), section 5.2:
  <https://wwwmpa.mpa-garching.mpg.de/~weiss/aa54956-25.pdf>.

Original Flaskamp and Braun sources were also read in the other checkout:
`/Users/owner/Documents/Software/dev/mesa_release/mesa_clean/main/notes/EbF/rsp2/references/`.
The six-page `flaskamp_numerical_scheme_note.pdf` is a later project summary,
not Flaskamp's thesis. Its equations were checked against the original.

The requested version retains kinetic-energy transport but uses LOCAL
relaxation closures for the other two moments. Setting their nonlocal
diffusion coefficients to zero without supplying these local closures would
give a different model. Braun's negative quantity was the temperature
gradient, not the temperature. Her solar experiments replaced one or both
closures only in the outer layers, recalibrated the models, and removed that
negative-gradient layer. They do not establish pulsation stability or show
that all nonlocal transport is unphysical.

## Definitions and two equations

Keep the current RSP2 specific turbulent energy:

```math
e_{\rm t}=w^2=\tfrac12\overline{v_i'v_i'}.
```

The new variables are the actual radial velocity--entropy covariance
`Pi = <v_r' s'>` and the full specific-entropy variance `<s'^2>`.
Flaskamp/Braun use `Phi = <s'^2>/2`; the equations below use the full variance
so its physical meaning and factors of two remain visible. `Pi` is NOT the
present code's `PII_face`: the latter has had a factor of velocity removed.

For uniform composition and the published isotropic local closure:

```math
\frac{D\Pi}{Dt}=
 \frac{2w^2}{3}\frac{c_p}{H_p}(\nabla-\nabla_{\rm ad})
 +\frac{T\nabla_{\rm ad}}{H_p}\overline{s'^2}
 -\left(\frac{\alpha_\Pi w}{\Lambda}+\frac1{\tau_{\rm rad}}\right)\Pi,
```

```math
\frac{D\overline{s'^2}}{Dt}=
 2\frac{c_p}{H_p}(\nabla-\nabla_{\rm ad})\Pi
 -\left(\frac{\alpha_\Phi w}{\Lambda}+\frac2{\tau_{\rm rad}}\right)
   \overline{s'^2}.
```

`D/Dt` follows the mean fluid. These are the reduced Kuhfuss moment equations,
not a claim to include all mean-strain and compressibility terms in a
pulsating star. No spatial transport of Pi or entropy variance appears.
Their local turnover damping is retained.

The luminosity and buoyant production use the same covariance:

```math
L_c=4\pi r^2(\rho T)_{\rm face}\Pi_{\rm face},\qquad
S_{\rm buoy}=\frac{T\nabla_{\rm ad}}{H_p}\Pi.
```

The second expression is continuous. A discrete cell source must be derived
from face moments and the EOS buoyancy coefficient, not copied from the old
`w_cell * average(PII/Hp)` while substituting a new Pi for PII.

Radiative cooling is

```math
\tau_{\rm rad}^{-1}=
 \frac{4\sigma\gamma_r^2 T^3}{c_p\kappa\rho^2\Lambda^2}.
```

`RSP2_alfar*x_GAMMAR` already supplies gamma_r, with
`x_GAMMAR = 2*sqrt(3)`. Zero alfar means zero inverse cooling time.
The current `compute_Dr` gives `w^2/tau_rad`. In the cited three-equation
model the cooling is in the entropy moments, with no extra direct Dr sink
in the turbulent-energy equation. Keeping both would be a different closure
and needs a derivation; it must not happen as an accidental carry-over.
The corresponding gas/turbulent exchange must be changed consistently in
both energy formulations and in LNA.

Braun et al. (2026), Table 1, list gamma_r=3.46, the rounded value of
2*sqrt(3). In RSP2 this is `RSP2_alfar=1d0`; the cooling rate scales as
`RSP2_alfar**2`. The one equation model evaluates cell Dr; the three equation
model evaluates face entropy cooling. Her alpha_omega=0.25 maps directly to
`RSP2_alfat=0.25d0` in `Lambda*sqrt(e_t)` kinetic energy diffusion, with no
extra sqrt(2/3) factor. Both existing RSP2 defaults remain zero. These are
coefficient mappings, not an assertion that the complete models coincide.

For any chosen temporal weight theta, the two residuals are

```math
R_\Pi=\Pi-\Pi_{\rm start}-\Delta t\left[
 (1-\theta_\Pi)(D\Pi/Dt)_{\rm start}+\theta_\Pi(D\Pi/Dt)\right]=0,
```

```math
R_{s^2}=\overline{s'^2}-\overline{s'^2}_{\rm start}-\Delta t\left[
 (1-\theta_{s^2})(D\overline{s'^2}/Dt)_{\rm start}
 +\theta_{s^2}(D\overline{s'^2}/Dt)\right]=0.
```

These formulas do not choose new time-centering controls. Time placement must
be specified with the existing pulsation scheme and reflected in LNA; making
the luminosity backward-Euler is not proposed.

## Parameters and viscosity

Two new local closure coefficients suffice for this restricted model:

| Coefficient | Term controlled | Published local MLT calibration |
| --- | --- | --- |
| alfa_pi | Entropy-flux relaxation | `6*sqrt(2/3) = 4.898979486` |
| alfa_phi | Entropy-variance relaxation | `4*sqrt(2/3) = 3.265986324` |

Implementation convention: these physical coefficients are source constants
`x_ALFAPI` and `x_ALFAPHI`, multiplied by `RSP2_alfa_pi` and `RSP2_alfa_phi`.
Both controls default to `1d0`. This recovers the local MLT calibration listed
in Braun et al. (2026), section 2; it is not her fitted solar calibration.
Her section 4.3 case A uses 2.88 and 2.0 when replacing both outer closures.
Her fiducial coefficients 2.155 and 2.0 multiply nonlocal transport;
rescaling those numbers does not calibrate our local decay controls.

Existing quantities provide Lambda, Hp, the kinetic dissipation coefficient
`C_D = RSP2_alfad*(8/3)*sqrt(2/3)`, gamma_r, kinetic-energy transport and
eddy viscosity. The old algebraic `x_ALFAS` source/flux prescription is
replaced by the new moment evolution in this mode, not retained as an extra
independent multiplier. No new nonlocal coefficients are needed for the
two added equations. A mode switch and complete state plumbing are separate
from the count of physical coefficients.

The published reduced equations omit explicit molecular-viscosity terms in
the two new moments. Their alfa terms are local closures of third-order
advection, and their radiative terms are thermal losses; they are not copies
of RSP2's eddy-viscous heating Eq. Existing Eq/Uq remain a paired exchange
between mean motion and turbulent kinetic energy.

This does not mean the exact entropy-flux equation has no dependence on mean
motion. Ahlborn's full RANS equation includes `-(Pi_vector dot grad) u`, whose
radial component for a radial mean flow is `-Pi * du/dr`. This is mean strain,
not an additional alpha_m dissipation term. The stellar-evolution reduction
omits it. Whether and how it is retained alongside RSP2's stress closure
must be resolved before calling a pulsation implementation complete. The
full equations also have correlations of viscous heating and fluctuating
pressure which the restricted three-equation model neglects.

## Y, EOS, and placement

The published driving term is `-ds/dr`. At uniform composition:

```math
-\frac{ds}{dr}=\frac{c_p}{H_p}(\nabla-\nabla_{\rm ad}).
```

Current RSP2 defines `gradT = gradL + Y_face`, and gradL may include
composition and dynamical modifications. In a coordinate where gradT is
the actual logarithmic temperature-to-pressure gradient, the difference is

```math
\nabla-\nabla_{\rm ad}
=Y_{\rm face}+\mathrm{gradL}-\nabla_{\rm ad}.
```

It is not automatically Y_face. The QHSE and radiation temperature rows
require their own coordinate mapping; substituting their gravity-normalized
gradT into this identity is not generally valid. In the exact entropy-gradient
expression, Hp is the pressure-gradient scale height. A dynamical implementation using
the existing gravity-based Hp helper must derive the relation consistently
from the selected pressure/temperature-gradient discretization. Simply
inserting a dynamical Ledoux threshold into the two published equations
does not accomplish this. Composition fluctuations require further physics
if a full Ledoux moment model is wanted; these two thermal moments alone do
not provide independent composition correlations.

The current EOS buoyancy factor is
`P*chiT/(rho*Cp*chiRho*Hp)`, equal to `T*grad_ad/Hp` for consistent
thermodynamics. Reuse this actual EOS expression and the public cell/face
length helpers when deriving the discrete equations.

Recommended placement after the user's staggering question: retain cell
e_t=w^2 and put BOTH new moments and their local evolution residuals at
faces. This is a discretization recommendation for the local version, not
a placement prescribed by the cited papers or an implemented choice.
Pi then supplies Lc directly at its existing location; entropy variance
couples directly to Pi at that same face. Putting variance in cells would
require interpolation in both directions between these two local equations.
Variance is an intensive second moment, not the cell's stored thermal
energy, so its name does not require cell placement.

At an interior face between cells k-1 and k, reconstruct the kinetic moment
from the cell energies with the selected RSP2 face weights:

```math
(e_t)_f=\mathrm{alfa}\,w_k^2+\mathrm{beta}\,w_{k-1}^2,
\qquad w_f^{\rm moment}=\sqrt{(e_t)_f}.
```

Use the same reconstructed moment in the face driving, local turnover
rates, and covariance bound. `average(w^2)` and `average(w)^2` differ on
an unequal state. The old luminosity uses `average(w)*PII`, which matters
for a flux-preserving conversion of an old snapshot. Regular derivatives
at vanishing turbulence still need explicit treatment; this formula alone
is not a solution to the zero-energy Jacobian issue.

The remaining face-to-cell coupling is buoyant production of cell turbulent
energy. Derive its average from the face covariance and EOS coefficient,
and use exactly the same cell exchange with the opposite sign in the gas
equation. Face Lc must enter the usual cell luminosity divergence. Conservation
depends on these identities, not on moving entropy variance into cells.
Boundary/forced-nonturbulent faces need explicit moment conditions compatible
with imposed Lc; they cannot blindly use every interior evolution row.
If nonlocal transport is later restored to the new moments, derive their
control-volume fluxes before retaining or changing this placement.

Y and L retain their existing gradient and flux-balance equations. Two new
moments mean two further unknowns and two further residuals, not replacing
Y or either luminosity/temperature-gradient row. Remesh, split/merge,
photos, retries, old state, boundary conditions, derivatives, and the two
additional LNA perturbations all need explicit handling.

## Analytic local equilibrium for startup

The following is a derived initializer for a quiet, homogeneous local
background: no kinetic-energy-flux divergence, turbulent pressure work,
mean strain or Eq forcing, and `gradT - grad_ad > 0`. It is not an exact
equilibrium of a pulsating or nonlocal envelope.

For general positive coefficients, define only the physically named rates

```math
\tau_\Pi^{-1}=\alpha_\Pi w/\Lambda+\tau_{\rm rad}^{-1},\qquad
\tau_{s^2}^{-1}=\alpha_\Phi w/\Lambda+2\tau_{\rm rad}^{-1}.
```

Holding w and the gradient fixed and setting only the two new time
derivatives to zero would give

```math
\Pi_0=\frac{(2/3)w^2(c_p/H_p)(\nabla-\nabla_{\rm ad})}
 {\tau_\Pi^{-1}
  -2(T\nabla_{\rm ad}/H_p)(c_p/H_p)(\nabla-\nabla_{\rm ad})
    /\tau_{s^2}^{-1}},\qquad
\overline{s'^2}_0=2\tau_{s^2}(c_p/H_p)(\nabla-\nabla_{\rm ad})\Pi_0.
```

This denominator can vanish or have the wrong sign. It is NOT an
unconditionally robust initializer for arbitrary frozen w and Y.

Instead solve all three local stationary equations together. At the
published coefficients, including `RSP2_alfad = 1`,
`2*alfa_pi = 3*C_D + alfa_phi`. The positive branch reduces to a quadratic:

```math
w_0^2+\frac{2\Lambda}{\alpha_\Phi\tau_{\rm rad}}w_0
=\frac3{16}c_pT\nabla_{\rm ad}
 \left(\frac{\Lambda}{H_p}\right)^2(\nabla-\nabla_{\rm ad}).
```

Use the rationalized positive root to avoid subtracting nearly equal large
numbers. If `w_ad^2` denotes the right-hand side (the no-radiation equilibrium
energy), then

```math
w_0=\frac{w_{\rm ad}^2}
 {\sqrt{[\Lambda/(\alpha_\Phi\tau_{\rm rad})]^2+w_{\rm ad}^2}
  +\Lambda/(\alpha_\Phi\tau_{\rm rad})}.
```

Recover the two moments from

```math
\Pi_0=\frac{C_Dw_0^3 H_p}{\Lambda T\nabla_{\rm ad}},\qquad
\overline{s'^2}_0=
 \frac{2(c_p/H_p)(\nabla-\nabla_{\rm ad})\Pi_0}
 {\alpha_\Phi w_0/\Lambda+2/\tau_{\rm rad}}.
```

When alfar=0 these simplify to

```math
w_0^2=\frac3{16}c_pT\nabla_{\rm ad}(\Lambda/H_p)^2
 (\nabla-\nabla_{\rm ad}),
```

```math
\Pi_0=\frac12\sqrt{\frac23}\Lambda w_0
 \frac{c_p}{H_p}(\nabla-\nabla_{\rm ad}),\qquad
\overline{s'^2}_0=\frac14\left[
 \Lambda\frac{c_p}{H_p}(\nabla-\nabla_{\rm ad})\right]^2.
```

Thus the local, no-radiation equilibrium matches the present algebraic
RSP2 luminosity closure in its corresponding homogeneous, hydrostatic
limit. Changing alfa coefficients or C_D generally removes the simple
quadratic reduction. The remaining local equation can be solved as a
scalar positive-root problem; negative or singular closed-form branches
must not be accepted. Stable zones with no imported turbulence use the
zero state; existing nonzero turbulence must not be erased merely because
the instantaneous local stratification is stable.

The prescribed gradient must also satisfy `L = Lr + Lc + Lt` for a static
envelope. Use these analytic moment relations inside that scalar gradient
solve and then the coupled static envelope solve. Independent local seeds
alone do not make the full envelope an equilibrium, particularly with
nonzero alfat or eddy-viscous production.

## Converting a pulsating snapshot

At an unforced face with positive entropy driving and positive old
convective luminosity, preserve that luminosity when initializing the
new covariance:

```math
\Pi_{0,f}=\frac{L_{c,f}^{\rm old}}{4\pi r_f^2(\rho T)_f}
          =w_f^{\rm old}\,PII_f^{\rm old}.
```

The missing entropy variance cannot be inferred uniquely from that snapshot.
A fully correlated parcel seed (consistent with the local MLT calibration)
is

```math
\overline{s'^2}_{0,f}=\frac{3\Pi_{0,f}^2}{2(e_t)_f},\qquad (e_t)_f>0.
```

Stable faces and faces without positive old luminosity start with
Pi=Phi=0, retaining their existing cell turbulent energy. This is a
one-time assumption for missing moments and changes the old heat flux at
those faces. It is not a reset of an evolved three equation state. The
target gradL must be evaluated and Y updated before classifying the face.
Otherwise a stale neutral gradient can manufacture entropy fluctuations.

If both energy and covariance vanish, initialize the variance to zero
without division. Nonzero covariance at zero energy is inconsistent and
must be handled in initialization, not hidden in a denominator floor.
This is a defined physical initial assumption, not a claim that the
instantaneous new residuals vanish or that the actual missing variance
has been recovered. An arbitrary old pulsation phase need not be a state
on the new model's periodic solution. Some physical adjustment is unavoidable.

The realizability condition is

```math
\Pi^2\leq\frac23 e_t\overline{s'^2},\qquad
e_t\geq0,\quad\overline{s'^2}\geq0.
```

The default local equilibrium above lies exactly on this covariance bound.
That is not proof that the evolution equations or a Newton solve preserve
it. For example, in the reduced isotropic equations without transport,
write `Q = (2/3)e_t<s'^2> - Pi^2`. At Q=0 and the published coefficients,

```math
\frac{DQ}{Dt}=-\frac43\overline{s'^2}\frac{De_t}{Dt}.
```

Thus positive local energy growth at that boundary can violate the assumed
isotropic covariance bound. This calculation is for the reduced equations,
not the complete Reynolds-stress system. It is a reason to audit
realizability/anisotropy before advertising this closure as a physically
enforced flux limit. Adding two equations is not by itself a guarantee of
positivity, realizability, or numerical stability. A variable transform or
clipping cannot repair an incompatible continuous closure.

## Implementation references and checks

- `star/private/hydro_rsp2.f90:compute_PII_from_Y`: algebraic PII being replaced
  only in a future three-equation mode.
- `compute_Lc_terms`: current face weights, `(rho*T)_face`, and `w_face*PII`.
- `compute_Source_div_w`, `compute_Source`: current factored source and EOS
  buoyancy expression; derive the new covariance-based cell source.
- `compute_D_div_w`, `compute_Dr_div_w`: existing C_D and inverse cooling time.
- `compute_RSP2_gradT`: actual gradL/Y relationship; dynamical and composition
  terms must not be silently equated with an entropy gradient.
- `do1_turbulent_energy_eqn`, `star/private/hydro_energy.f90`: energy exchange,
  work, and time-centering consistency.
- `star/private/star_LNA_turbulence_closures.f90` and
  `star/private/star_LNA_support.f90`: future independent moment perturbations;
  the current algebraic perturbation is insufficient for a three-equation mode.

Analytic verification on 2026-09-19 used standalone Python, not MESA:
108 combinations of buoyancy coefficient, entropy gradient, mixing length,
and radiative cooling rate satisfied all three local stationary equations;
maximum residual divided by the sum of term magnitudes was 3.53e-16.
The same cases confirmed saturation of the isotropic covariance bound.
These checks verify the initializer algebra, not stellar-model behavior.

Outstanding before implementation: select the pulsation mean-strain
approximation; derive entropy-gradient/buoyancy discretization with the
existing dynamical options; decide face reconstruction and source averaging;
verify continuous realizability; specify time centering and full state/LNA
integration. This proposal does not resolve the separate existing k+2
conservative-work Jacobian issue.


## Discretization follow-up, 2026-09-19

Sections 2.1--2.5 of the [implementation plan](rsp2_three_equation_implementation.md)
record the temperature-row entropy-gradient audit, consistent face EOS
coefficients, luminosity-based moment scales, and a cell-source reconstruction
using face Pi/sqrt(e_t_face) with the cell's own w. The initial logarithmic
entropy-gradient inversion has been withdrawn: a follow-up check shows
finite-zone bias at Y=0 even with exact discrete hydrostatic balance. It
can dominate a small evolutionary driving term. The replacement is now
implemented in the inactive three-equation path: it uses each temperature
row's discrete thermodynamic differential and sets its gradL reference and
thermal driving together. The actual-log row needs no HSE-coordinate
conversion; standard and dPrad use their actual row coefficients. Composition
offsets and the actual-versus-HSE pressure term remain explicit. Standalone
neutral-state, tiny-Y, reference-invariance and derivative checks pass;
actual EOS/Fortran AD and MESA runs remain untested. See section 2.1 of the
implementation plan for equations, source references and precision limits.
Separately, the cell-source reconstruction replaces the provisional direct
average of face buoyancy products,
which can deplete a zero-energy cell adjacent to finite face turbulence.
The source reconstruction is a numerical choice consistent with the
homogeneous continuous limit, not a new flux cap or a prescription quoted
from these papers. The plan records the local canceled energy residual,
its limits for nonzero Lt/u-flag face heating, and the tested two-cell
prototype. Existing Eq/Uq and Lt formulas are retained. These standalone
checks do not establish full MESA convergence or continuous realizability.
