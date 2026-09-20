# RSP3 covariance closure correction

Date: 2026-09-20. Status: derivation and standalone checks, not implemented.
The two requested stellar restarts are complete. No production source was
edited; all test input edits were confined to isolated copies. The user authorized deriving
a consistent correction after the variance failure was identified.

The subsequent [placement audit](rsp3_layout_audit.md) compares the current
staggering, cell moments and face moments. Its preferred simple long-term
RSP3 design places all three moments on faces, after rederiving the energy
projection. It identifies a thermal coupling defect in the simple cell
averaging candidate below. The face layout has since been implemented and
installed; the local covariance correction below remains unimplemented.
The local closure derivation remains applicable to moments evaluated together.
The subsequent saved-output audit is in [rsp3_face_profile_audit.md](rsp3_face_profile_audit.md).
The current finite-step derivation, local/nonlocal design and validation record
are in [rsp3_covariance_discretization.md](rsp3_covariance_discretization.md).
The user has clarified that both local and nonlocal operation are required,
including long evolutionary timesteps and moving convective boundaries.

## Current recommendation

Retain the common face placement of w, Pi and Phi. It keeps Pi at the
enthalpy-flux interface and lets local sources use colocated moments.
The placement audit demonstrates the lost alternating thermal response
of the simple all-cell averaging alternative. It does not prove that
faces outperform every possible cell formulation. The present local
closure failure also exists without any mesh.

Proceed with the native implementation design specified in
[rsp3_covariance_discretization.md](rsp3_covariance_discretization.md).
Treat the one-third buoyancy proposal as a
candidate, not as a recovered literature coefficient or a verified fix
for every stellar failure. Establish the pressure/entropy approximation,
decay and compression together, then verify the complete residual and AD
derivatives, initialization, work and transport. The agreed scope includes
both local and nonlocal operation. Nonzero alfat uses common transport of
all three moments. Preserve the current radiative and convective heat-flux
time weighting; the finite-step proposal treats moment transport implicitly
and matches that Lt time level in gas energy. Enhanced stable-layer
dissipation is deferred.

Native comparisons must distinguish moment-domain preservation, residual
convergence, total energy conservation and pulsation growth. Passing the
first property alone does not establish the other three. Update LNA and
the manuscript from the accepted equations after that derivation.

## Run evidence

### Face layout in star/work2, 2026-09-20

The user's later 15 Msun pre-main-sequence run provides an independent
recurrence after the face conversion. The supplied terminal trace ends
after 576 accepted models, with 111 retries. The largest final residual
on rejected attempts is dPhi_dt 34 times, dPi_dt 69 times, detrb_d 6 times
and rsp2_fl twice. The run advances to about 925.85 years; it is not
permanently stalled, but repeatedly spends 25 iterations on fixed residuals.
Median accepted iterations are five; the 90th percentile is ten.

The current case uses v_flag, dedt with conservative work, dynamical gradL,
alfa_pi=alfa_phi=1, and alfat=alfam=0. The defaults leave alfar=0.
The face conversion did not apply the local covariance correction below.

At model 56, solver call 391, dPhi_dt at face 451 stalls at 1.3278e-3
through iteration 25. At models 100, 150 and 200, the same face is accepted
at iteration 10 only after the maximum residual tolerance becomes 1e-4.
Its residuals remain 5.6777e-5, 5.4364e-5 and 5.9869e-5 respectively.

The saved version-22 photos independently show:

| Model | Face | Y | w, cm/s | Pi | Phi |
| --- | --- | --- | --- | --- | --- |
| 50 | 451 | 1.79026e-4 | 1.03633e4 | 6.46942e7 | 5.19206e7 |
| 100 | 451 | -2.19960e-4 | 1.17728e4 | 2.40284e7 | 0 |
| 150 | 451 | -1.94475e-4 | 1.18051e4 | 1.80749e7 | 0 |
| 200 | 451 | -1.31935e-4 | 1.18203e4 | 8.47659e6 | 0 |

Pi and Phi have the units defined below. Nonzero Pi with zero full entropy
variance is inadmissible regardless of the isotropic velocity assumption.
Here composition is initially uniform and dynamical gradL is enabled,
so `get_rsp2_thermal_gradient` gives the sign of -ds/dr from the signed Y
times its positive temperature-gradient coefficient. The negative Y and
positive Pi therefore drive Phi further negative when Phi is already zero.
`solver_support:Bdomain` clips the negative Phi update. Additional Newton
iterations cannot repair the underlying source/domain conflict at this state.

Holding the other coefficients at a trial state, the variance row asks for

```math
\Phi=\frac{\Phi_{\rm start}-2\Delta t\Pi\,\partial s/\partial r}
 {1+\Delta t[C_\Phi w/\Lambda+2/\tau_{\rm rad}]}.
```

A negative numerator gives an inadmissible trial variance. Changing a
positive residual normalization cannot change that sign. Additional decay
proportional to Phi, including radiative cooling or a buoyancy decay rate,
vanishes at Phi=0 and does not by itself fix this boundary failure.

Later stalled dPi_dt attempts are a separate unresolved part of this run.
At model 575, face 428 stays at 4.9254e-3 while corrections are approximately
1e-11; the smaller-timestep retry converges in five iterations. At model
576 the corresponding plateau is 3.3876e-3 at face 444. Saved accepted
photos do not contain those failed trial states. Do not assert that all
69 dPi_dt retries are proved to have the same cause as the early Phi failure.

Reproducible parser, copied inputs/photos, full attempt records and summary:
`output/review/rsp3_work2_moments_20260920`. Photo record markers, variable
indices, model numbers and layout version are checked. The read-only audit
did not compile or run MESA, change the case, change tolerances or implement
a closure correction.

### Earlier staggered layout restarts

Both runs restart the user's executable and photo x00115000 in isolated
copies, stop at model 116500, and save every accepted profile. The mesh is
fixed at 149 zones. The only physics difference between them is alfat.

| alfat | Accepted steps | New retries | Minimum dt, s | Median accepted iterations | Mean accepted iterations |
| --- | ---: | ---: | ---: | ---: | ---: |
| 0 | 1500 | 80 | 12.6542 | 4 | 6.3033 |
| 0.1 | 1500 | 87 | 11.3893 | 4 | 6.1747 |

Every retry has dPhi_dt as the largest residual, at face 74 or 75.
Accepted constrained steps can also stall at face 76. At model 115550,
face 76 has Phi=0 and Pi=1.43407e11. Reconstructing the variance residual
from the saved profile gives 3.81799136e-5, matching the printed 3.8180e-5.
The residual does not fall during iterations 3 through 10. It is accepted
when the existing tolerance changes to 1e-4. Increasing dt eventually
pushes this same residual above that tolerance and causes another retry.

The relevant physical row is

```math
R_\Phi=\Phi-\Phi_{\rm start}-\Delta t\left[
 2\left(-\frac{\partial s}{\partial r}\right)\Pi
 -\left(\frac{C_\Phi w}{\Lambda}+\frac{2}{\tau_{\rm rad}}\right)\Phi\right].
```

At Phi_start=Phi=0, positive Pi and a negative entropy driving gradient,
this residual is strictly positive for every positive timestep. The
Newton update asks for negative variance and Bdomain clips that update.
Changing residual normalization cannot make the physical row vanish.

The full run outputs, executable/photo hashes, JSON statistics, profile
reconstruction, local ODE checks and figure are in
`output/review/rsp3_restart_115000_20260920`. `analyze.py` reproduces the
analysis. These runs do not establish how a corrected model will behave.

### Face layout in the pulsation case, 2026-09-20

The later pasted trace covers models 6319--6463. Of these 145 attempts,
73 take five iterations and 72 take ten. Models 6336--6407 all stall on
dPhi_dt at face 198. For example, model 6337 stays at 1.5920e-5, and model
6372 at 6.8780e-6, from iteration two through ten. They are accepted when
the maximum tolerance becomes 1e-4. There is no new retry in this excerpt;
the printed cumulative count remains 14.

The exact failed/intermediate trial states are not saved. However, the
subsequent photo at model 7000 independently demonstrates the same domain
failure at that face: w=1.32363e5 cm/s, Y=-0.0226495, Pi=3.92819e9, Phi=0,
and Lc/L=5.50423e-5. The photos at 5000 and 8000 have covariance ratios
Pi^2/[(2/3)w^2 Phi] of 1.3138 and 1.2553 at face 198. Thus the pulsation
model also leaves the moment domain; the trace alone was not used to infer
the value of Phi at model 6337. Its current run differs from the earlier
model-1000 profile audit, and copied inputs are kept separate.

The parser and copied photos are in
`output/review/rsp3_work2_moments_20260920/pulsation`; the parent directory's
`analyze_pulsation.py` reproduces the summary without MESA.

## What must remain possible

Pi is a signed velocity/entropy covariance. Negative Pi is allowed.
Positive Pi in a layer with a negative entropy driving gradient is also
allowed. It describes a countergradient heat flux.

Phi is the full entropy variance, `<s'^2>`, and must be nonnegative.
Zero Phi requires zero Pi. Under the adopted isotropic velocity closure,

```math
\Pi^2\leq\frac23e_t\Phi,\qquad e_t=w^2.
```

This is a condition on actual second moments, not an empirical cap on
convective luminosity. Changing the representation of Phi cannot make
an evolution law that leaves this domain physically consistent.

## Local source correction

### Published coefficients and the implemented normalization

A direct check against Braun et al. (2026), section 2, equations (2)--(4)
and (11), confirms that the implemented isotropy factors were not omitted.
The same factors appear in Kupka, Ahlborn and Weiss (2022), Appendix A.
Their kinetic energy omega is our w^2. Their Phi is half our full variance.
Consequently their `2*(grad_ad*T/Hp)*Phi` becomes
`(grad_ad*T/Hp)*Phi` in our Pi equation, and their variance source doubles.
The implemented entropy-gradient source in Pi is already
`(2/3)*w^2*(-ds/dr)`. The local decay constants `6*sqrt(2/3)` and
`4*sqrt(2/3)` also match the published local MLT calibration.

The factor one third proposed below is therefore a new closure choice.
It is not a recovered isotropy factor, a correction for the full-variance
normalization, or a demonstrated transcription error in the published
buoyancy term. Isotropy alone does not select this modified pressure/entropy
closure. The moment-domain counterexample motivates investigating such a
closure, but does not establish that this particular choice is the best
physical model for pulsation. The earlier statement that the production
term simply needed to be made isotropic was too categorical.

### Proposed modified closure

For this derivation, all three moments are evaluated at the same location.
The buoyancy coefficient is the existing EOS expression
`-(dP/dr)*chiT/(rho*chiRho*Cp)`. The entropy driving is the existing mapped
`-ds/dr`, including the selected temperature row. Neither is replaced by
an HSE-only approximation.

Keep the turbulent energy equation's local source and dissipation:

```math
\frac{D e_t}{Dt}=
 -\frac{1}{\rho}\frac{\partial P}{\partial r}
   \frac{\chi_T}{\chi_\rho c_p}\Pi
 -\frac{C_D}{\Lambda}e_t^{3/2}+E_q
 -\frac23\alpha_p e_t\,\nabla\!\cdot u.
```

Here C_D=RSP2_alfad*(8/3)*sqrt(2/3), and
C_Phi=RSP2_alfa_phi*4*sqrt(2/3). The variance equation remains

```math
\frac{D\Phi}{Dt}=
 2\left(-\frac{\partial s}{\partial r}\right)\Pi
 -\left(\frac{C_\Phi w}{\Lambda}+\frac{2}{\tau_{\rm rad}}\right)\Phi.
```

A covariance-consistent isotropic flux equation is

```math
\frac{D\Pi}{Dt}=
 \frac23e_t\left(-\frac{\partial s}{\partial r}\right)
 -\frac{1}{3\rho}\frac{\partial P}{\partial r}
       \frac{\chi_T}{\chi_\rho c_p}\Phi
 -\left[\frac{C_D+C_\Phi}{2}\frac{w}{\Lambda}
        +\frac{1}{\tau_{\rm rad}}
        +\frac{\alpha_p}{3}\nabla\!\cdot u\right]\Pi.
```

There are three linked changes, not just a tuned decay coefficient:

1. The flux buoyancy term has one third of the old coefficient. The
   radial velocity variance is 2e_t/3, so its buoyancy production is
   2/3 of the energy production. The paired flux production must use
   the compatible coefficient. This represents isotropic redistribution
   of buoyancy production and modifies the reduced Kuhfuss closure.
2. The flux decay is half the sum of the kinetic-energy and entropy-
   variance decay rates. With unit controls it is (10/3)*sqrt(2/3),
   replacing 6*sqrt(2/3).
3. Mean deformation must be treated with the same isotropic compression
   as the energy equation. Keeping `-Pi*du/dr` alone while imposing an
   isotropic kinetic moment, or disabling pressure work with alfap=0,
   does not give the same covariance evolution. A model retaining the
   distinct radial strain needs a corresponding radial stress closure.

For these local terms, direct differentiation gives

```math
\frac{D}{Dt}\left(\frac23e_t\Phi-\Pi^2\right)=
 -\left[\frac{(C_D+C_\Phi)w}{\Lambda}
       +\frac{2}{\tau_{\rm rad}}+\frac23\alpha_p\nabla\!\cdot u\right]
      \left(\frac23e_t\Phi-\Pi^2\right)
 +\frac23 E_q\Phi.
```

Thus the local covariance boundary cannot be crossed outward if E_q is
nonnegative. The implemented viscous heating and any time discretization
of it must be checked against that assumption; it is not an assumption
about arbitrary user energy sources.

Additional nonnegative flux decorrelation adds a nonnegative term to the
right-hand side. One possible interpretation of RSP2_alfa_pi is to
multiply the common turnover decay `(C_D+C_Phi)/2`, with alfa_pi>=1.
Its default remains one. This changes the control's meaning and would
require explicit documentation and validation of coefficient combinations.
Do not retain the old independent multipliers without checking the bound.

This is one internally consistent isotropic closure. The covariance
identity establishes consistency, not a calibration of turbulent pressure
correlations or a unique physical closure. It is not claimed to reproduce
the original model's transient growth rates.

## Local equilibrium and continuous checks

For uniform composition, no radiation, mean strain, transport or viscous
heating, and unit controls, both the old and corrected equations have the
same positive stationary solution. In dimensionless units with buoyancy,
positive entropy driving and mixing length equal to one,

```math
e_t=0.1875,\qquad\Pi=0.176776695296637,\qquad\Phi=0.25.
```

Equivalently the dimensional energy is

```math
e_t=\frac{3}{16}\Lambda^2
 \left[-\frac{1}{\rho}\frac{\partial P}{\partial r}
 \frac{\chi_T}{\chi_\rho c_p}\right]
 \left(-\frac{\partial s}{\partial r}\right).
```

The source reduction and decay change must be made together for this
agreement. It is an equilibrium check, not agreement of the time-dependent
solutions.

An initially admissible homogeneous stable state uses e_t=1e-12, Pi=0,
Phi=1, buoyancy=1, entropy driving=-1 and Lambda=1. The old unit-coefficient
ODE crosses Phi=0 at dimensionless time 1.4165846007 with Pi=0.0631542
and dPhi/dt=-0.1263084. No mesh, Newton method, radiative loss, transport
or mean deformation is present. The corrected ODE stays nonnegative
through time 10 and retains a nonnegative covariance determinant.
Both ODE tolerances, 1e-9 and 1e-11, agree. A further 1000 random states
verify the determinant identity, including compression, radiative cooling,
nonnegative viscous heating and extra decorrelation; maximum scaled error
is 4.76e-14. These are local checks only.

## Implicit step and implementation contract

This section specifies the numerical work required in addition to the
continuous closure. It is a proposal, not an implemented or stellar-tested
fix. The source, gas energy equation and case controls remain unchanged.

The covariance matrix is a mathematical representation of the existing
three moments, not an additional solver variable:

```math
\begin{pmatrix}(2/3)w^2&\Pi\\\Pi&\Phi\end{pmatrix}\succeq0.
```

The corrected local sources evolve this matrix by a common linear drift
at fixed trial gas state and turnover rates, plus a nonnegative variance
source from Eq. The off-diagonal drift coefficients are the existing
`-(dP/dr)*chiT/(3*rho*chiRho*Cp)` and the existing mapped `-ds/dr`.
The diagonal amplitude decay rates are half the respective variance
decay rates. This is why the same production and loss cannot be chosen
independently in the three equations.

For stable stratification, no expansive drift from compression and
nonnegative Eq, the fixed-coefficient drift is stable. Its backward Euler
resolvent is an integral of positive covariance evolutions, so it preserves
the domain for every positive timestep. Nonlinear turnover damping can
still be evaluated implicitly. This does not require making L fully
implicit or changing its theta weighting.

The standalone `check_closure.py` in
`output/review/rsp3_work2_moments_20260920` now checks:

- The old continuous ODE crosses Phi=0 at time 1.41658460065 with Pi=0.0631542;
  the corrected ODE remains admissible through time 10.
- 3000 stable backward Euler blocks, timesteps from 1e-6 to 1e6 and
  nonnegative driving and decorrelation. The smallest relative covariance
  eigenvalue is 4.96e-8; maximum scaled linear-equation error is 2.40e-16.
- 300 nonlinear stable backward Euler problems with the rates proportional
  to the new w. The energy consistency error is at most 1.07e-11 relative.
- A fixed-coefficient unstable example in which backward Euler at too large
  a timestep gives negative variances. Continuous realizability alone is
  therefore not an unconditional guarantee for every flow and timestep.
- A two-zone energy-only diffusion counterexample and the admissible result
  obtained with one common implicit diffusion operator.
- The unchanged unit-control, uniform convective stationary state.

The local scalar root calculation used in the nonlinear check is a test,
not a proposal to replace the global MESA solver. These tests establish
the stated local properties, not global Newton convergence in a star.

The actual `do1_turbulent_energy_eqn` uses `w^2-w_start^2`.
Its divided row is the same energy equation divided by the new positive w;
it is not backward Euler on w itself. Thus the local backward Euler checks
use the correct conserved moment. `calc_Ptrb_work_face` additionally contains
both new and old energy when pressure time centering is enabled. The
continuous compression formula alone does not cover that discrete work.

There is also a startup predictor inconsistency. In
`RSP2_adjust_vars_before_call_solver`, starting from w=Pi=0 and Phi>0,
the current predictor sets w=dt*abs(buoyancy)*sqrt(Phi) and
Pi=dt*buoyancy*Phi. Its covariance ratio is exactly 3/2, above the allowed
value one. This is a trial-state defect, not evidence that this branch
caused the later face-198 plateau. It must be corrected along with the
source. A coupled moment predictor can remain admissible; changing only
the source coefficients while keeping this predictor is incomplete.

The implementation contract is:

1. Keep w, Pi and Phi on faces. Change the Pi buoyancy coefficient, turnover
   decay and deformation consistently with the energy and variance rows.
   The default alfa_pi remains one, but its documented normalization becomes
   a multiplier of the minimum compatible decay. Enforce alfa_pi>=1 for
   this interpretation, with nonnegative alfad and alfa_phi. This is a new
   isotropic closure, not a claim that the old physical model had a missing
   factor in its Fortran transcription.
2. Retain the implicit moment time level and the current heat-flux formula.
   Preserve all gas-state and w derivatives. Check the complete three-row
   Newton block against finite differences. An algebraically equivalent
   local block preconditioner is available if stiffness remains, but changing
   residual units alone cannot correct an inadmissible physical source.
3. Use admissible predictors and check the full covariance domain in Newton
   globalization and before accepting a model. Independent clipping of Phi
   is not a complete domain treatment. The exact-zero state and rank-one
   covariance boundary need explicit tests; do not introduce a numerical
   variance floor. A rejected physical branch must trigger a changed trial
   or retry, not become accepted solely because tol3 is looser. Do not cap
   all evolutionary steps at the turnover time: use the nonlinear saturated
   branch where it exists, with admissibility checked at the trial state.
4. For turbulent pressure, use the same fractional energy compression as
   the actual discrete w row, with half that rate in the Pi amplitude.
   The continuum alpha_p*div(u)/3 expression is its limiting form. Reusing
   an unrelated radial strain would invalidate the paired source proof.
   Check both u/v and both gas-energy work forms with their existing
   cell/face projection; do not change total energy accounting independently.
5. Nonnegative viscous energy input is part of the covariance proof.
   Current Eq uses new strain times time-centered strain and can be negative
   during a reversal. A general dissipative design must use the same strain
   in the stress and its work, updating Uq and Eq together so that work is
   nonnegative and the discrete mechanical-energy identity still holds.
   This temporal choice requires validation; simply clipping negative Eq
   would violate energy conservation. TDC must remain separate. This issue
   cannot explain work2's alfam=0 failure.
6. For nonzero alfat, use the common moment transport closure described
   below. Leaving energy-only diffusion enabled while claiming general
   covariance preservation is not a complete solution. The user has since
   clarified that both local and nonlocal operation are required; the
   finite-step proposal specifies their shared operator and time level.
7. The common positive overlap remap already preserves admissibility when
   the donor moments are admissible. Test ordinary remesh, split/merge and
   envelope remesh after the source change. Invalid old photos must not
   silently pass as valid corrected initial states or be silently repaired;
   prepare and relax a consistent background explicitly.
8. Update star_LNA from these same physical source derivatives. The local
   stable oscillation frequency changes from approximately sqrt(8/3)*N to
   sqrt(4/3)*N before damping. Physical moment modes need not disappear;
   mode selection and eigenvector accuracy remain separate checks. Validate
   the radial fundamental by shape and convergence, not table index alone.

This is a complete set of design obligations, not proof that every listed
numerical choice has passed a stellar test. In particular, compression,
viscous time weighting, full Newton globalization and nonlocal transport
cannot be certified by the local ODE checks. Compare both supplied cases
at fixed physics and successively smaller timesteps, verify the covariance
bound and actual residuals, and measure energy error and pulsation growth
before accepting the change. A physical closure choice precedes coding.

## Nonzero alfat

Transporting only energy does not preserve the covariance domain. If a
face or cell is on its covariance boundary and energy is exported,

```math
\left.\frac{D}{Dt}\left(\frac23e_t\Phi-\Pi^2\right)\right|_{L_t}
 =-\frac23\Phi\frac{\partial L_t}{\partial m},
```

which is negative for a net energy loss. The other two moments cannot
remain unchanged under that export while retaining the same covariance.
This is independent of how small or positive alfat is.

A sufficient consistent transport closure is to apply the same positive
mass diffusion operator to all three moments:

```math
\left.\frac{D f}{Dt}\right|_{\rm transport}
 =\frac{\partial}{\partial m}\left[
 (4\pi r^2\rho)^2\,\alpha_t\Lambda w\,
 \frac{\partial f}{\partial m}\right],
 \qquad f=e_t,\Pi,\Phi.
```

The energy part is the existing Lt form. The coefficient and face weights
must be identical in all three rows. A common monotone finite-volume
operator mixes admissible covariance matrices with nonnegative weights.
Different independent diffusion coefficients do not give this guarantee.
This is sufficient, not the only possible realizable transport closure.

This proposal adds Pi and Phi transport. The earlier local-only scope was
superseded by the user's request to support both local and nonlocal models.
Common moment transport does not guarantee the absence of a negative mean temperature
gradient; realizability and thermal stratification are different tests.

## Cell and face placement

The earlier recommendation to put Pi/Phi on faces prioritized the direct
Lc coupling but did not establish compatibility with cell kinetic energy.
It should not be read as a demonstrated superiority of face placement.
Likewise, the cell proposal below is a candidate, not a demonstrated
convergence improvement or an approved layout change. The homogeneous
closure failure exists without any spatial placement. Establish the local
closure first, then compare compatible cell and staggered discretizations.
The later placement audit supplies this comparison and supersedes a
preference for the simple cell reconstruction described below.

The local proof does not apply directly to today's cell w plus face Pi/Phi.
The face energy is an average of two cell energies, and each cell source
averages its two bounding Pi values. Consequently its reconstructed face
energy derivative includes neighboring entropy fluxes, while its Pi row
uses only its own Phi. Even with equal cell energies and weights, setting
both neighboring fluxes to zero makes the central face's energy production
half the local expression used in the covariance derivation. On its
covariance boundary, the corrected local buoyancy coefficients would then
still give a negative determinant derivative for positive central flux.
Changing only the two constants in hydro_rsp2 is therefore insufficient.

A two-cell check gives a covariance determinant of +0.014 after common
positive mixing, compared with -0.18 after the identical energy transfer
without Pi/Phi transport. The equal-weight staggered buoyancy example
gives determinant derivative -0.272165527. These explicit counterexamples
are in `transport_placement_checks.json` alongside the run analysis.

A simple candidate for a complete formulation keeps the existing cell e_t and stores
both Pi and Phi in those same cells. Face Lc still uses
`4*pi*r**2*rho_face*T_face*Pi_face`, with a common convex reconstruction
of e_t, Pi and Phi to the face. An average of admissible covariance
matrices is admissible. The independent face Y and luminosity equations
retain their locations. Pi_face would be a derived value, not a new solver
variable. Current face moment photos cannot silently be reinterpreted as
cell data, and this failing photo already contains incompatible moments.
A documented conversion/initialization is required before a corrected run.

### Face luminosity residual with cell moments

Cell storage does not remove face fluxes or the independent Y equation.
For an interior face k between cells k-1 and k, the simplest common convex
reconstruction is

```math
\Pi_{f,k}=\mathrm{alfa}_k\Pi_k+(1-\mathrm{alfa}_k)\Pi_{k-1},
```

and the existing face residual remains

```math
R_{{\rm flux},k}=\frac{
 L_{r,k}(Y_{f,k})+4\pi r_k^2\rho_{f,k}T_{f,k}\Pi_{f,k}
 +L_{t,k}-L_k}{L_{{\rm scale},k}}=0.
```

All current thermodynamic and reconstruction dependencies must remain in
the Jacobian. With these fixed for the partial derivative, the two cell
Pi columns have coefficients `area*rho_face*T_face/L_scale` times their
respective reconstruction weights. The ordinary diffusive radiative row
retains `dR_flux/dY_face=Lrad_coeff/L_scale`. Thus the face row continues
to constrain Y, and each cell Pi is independently constrained by its
moment evolution row. No additional independent face Pi is required.
The surface retains its separate existing luminosity boundary condition.

The same computed Lc at an internal face must enter both adjacent energy
balances with opposite signs and the existing temporal weights. The
additional nonlocal transport of Pi and Phi also uses shared face fluxes
computed from the neighboring cell moments. It is distinct from Lc.

This establishes equation coupling and conservation, not Newton convergence.
On a uniform mesh, averaging an exactly alternating cell Pi perturbation
gives zero interior face Pi. Its moment storage and decay rows still
constrain that perturbation, and common diffusion damps it when enabled.
Check the full coupled Jacobian and the alfat=0 limit before claiming that
this reconstruction is sufficiently accurate or well conditioned for the
stellar run. Preserving local covariance does not prove this spatial or
temporal discretization is the best one. A compatible staggered formulation
remains an alternative if it can be derived and verified.

The subsequent audit demonstrates that averaging face entropy gradients
into cells and Pi back to faces loses the response to alternating mean
entropy in a constant-coefficient test. Pi storage and decay do not repair
that defect, because this entropy perturbation produces no Pi perturbation.
Thus this simple reconstruction is not recommended as a complete scheme.

For remeshing, use the same positive overlap weights for all three cell
moments; conserve cell turbulent energy and preserve admissibility. The
current separate face interpolation of Pi/Phi is not appropriate for this
candidate. Higher-order independent scalar reconstructions are not by
themselves covariance preserving.

## Required implementation and verification

For the cell candidate if a suitable thermal flux discretization is derived
and that formulation is selected; the later audit records separate
requirements for the preferred all-face candidate:

- `hydro_rsp2:rsp2_moment_rhs`, `rsp2_moment_source` and
  `do1_rsp2_moment_eqns`: colocated production, paired decay/compression,
  and common transport with complete AD derivatives.
- `compute_Lc_terms`: reconstruct Pi to its existing face luminosity.
  Preserve the enthalpy and total-energy flux accounting.
- `solver_support:Bdomain`: maintain an admissible trial state; do not
  accept a clipped state merely because the residual scale is large.
  Derive the discrete time update as well as the continuous identity.
- State/photo/model layout, initialization, three remesh paths and profile
  column locations: explicit conversion and the same covariance convention.
- LNA: linearize the actual corrected closure and transport; update work
  terms and both manuscripts. The equilibrium match does not imply an
  unchanged eigenvalue growth rate.
- Validate local decay/growth, both gradient signs, zero energy/variance,
  nonzero alfat, pressure work, u/v heating, all temperature rows, restart,
  ordinary mesh adjustment and split/merge AMR. Then repeat the two full
  stellar runs and assess growth rates and conservation separately.

The source equations behind the original local calibration are Braun et al.
(2026), section 2 and section 4.3,
https://arxiv.org/html/2604.06151v1 . The proposed covariance correction and
its consistency checks above are derived here; they are not attributed to
that paper.
