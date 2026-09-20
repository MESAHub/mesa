# RSP3 moment placement audit

2026-09-20. Design audit and standalone checks only. No production source,
inlists or saved models changed. No MESA compilation or new stellar run.

The comprehensive conversion plan is now in
[rsp2_face_w_implementation.md](rsp2_face_w_implementation.md). It covers
the residuals, gas energy projection, TDC-style eddy-viscosity discretization
in separate routines, all mesh paths, saved state, LNA and validation. The
eventual target is common face energy machinery for RSP2 and RSP3, after
separate validation; the evidence below first motivates RSP3 collocation.

## Recommendation and limits

For a three-equation model intended for pulsation and stellar evolution,
including later nonlocal moment transport, the preferred simple design from
this audit places all three moments on faces: e_t=w^2, Pi and Phi. The gas
thermodynamic state remains in cells. Y and total luminosity remain on faces.
This is a recommendation for the eventual RSP3 formulation, not an approved
or implemented change to RSP2's existing cell w.

The reasons are coupled, not a general rule that faces are better:

1. Pi remains where the temperature difference drives heat transport and
   Lc enters the two neighboring energy balances.
2. All three moments share the location at which their covariance condition,
   sources, decay and future nonlocal transport must be consistent.
3. Common conservative diffusion and remapping can operate on the same
   face control volumes for all three moments.

The current mixed layout does not supply the second property. The simple
cell layout previously proposed does not supply the first adequately:
averaging gradients into cells and Pi back to faces leaves an alternating
temperature pattern undetected by convective heat transport. A more
sophisticated cell formulation can avoid that defect. It has not been
derived or compared here, so this audit does not establish that face
moments are superior to every possible cell method.

All-face moments require more changes to RSP2 energy accounting than cell
moments. This cost is material. The recommendation concerns the eventual
coupled equations, not the smallest diff. It does not establish fewer
Newton iterations in MESA. That requires the verification below.

## Correction to the previous advice

The original recommendation put Pi and Phi on faces because Lc is a face
flux. It omitted the consistency check with cell e_t. The subsequent cell
proposal prioritized local covariance, transport and remeshing. It omitted
the steady thermal response after interpolating both gradients and fluxes.
Both recommendations were incomplete.

The cell Pi storage and decay row does constrain an alternating Pi
perturbation. But an alternating *mean entropy* perturbation can produce
no cell gradient and therefore no Pi perturbation. The moment row then
has nothing to damp. Adding Pi diffusion does not repair that absent
thermal coupling. The full MESA Jacobian need not be singular: radiation,
coefficient variations and other equations contribute. The demonstrated
defect is a missing convective response.

No placement change cures the homogeneous closure failure in
[rsp3_covariance_closure.md](rsp3_covariance_closure.md). That ODE has no
mesh. The proposed isotropic correction remains a physical model choice,
not a missing factor established from the original literature.

## Stable-layer face flux audit, model 1001

The current 350-zone RSP3 snapshot gives direct evidence of the averaging
weakness. `hydro_rsp2:rsp2_moment_source` contributes half of each bounding
face's buoyancy times Pi, multiplied by cell w divided by that face's RMS w.
For uniform background coefficients and w this reduces to

```math
 \left.\dot e_{t,k}\right|_{\rm buoyancy}
 \mathrel{\propto}\tfrac12(\Pi_k+\Pi_{k+1}).
```

An alternating face Pi perturbation produces zero cell energy forcing.
Reconstructing face energy from those cell energies averages a second time.
On a uniform mesh its response is multiplied by cos(theta/2)^2. It vanishes
at the alternating-grid pattern. The face Pi/Phi rows and the gas heating
from the difference of adjacent Lc values still respond to that pattern.
This is a loss of the local kinetic-energy response, not a claim that the
full stellar Jacobian is singular. End conditions and variable coefficients
also prevent treating the periodic example as an exact global null vector.

In the saved model, between 1 and 3 Rsun there are 24 Pi sign changes across
37 stable faces. Six cell source sums cancel by more than 80 percent. At
cell 295, r=1.25347 Rsun, the two buoyancy contributions are -1.647106e-4
and +1.647768e-4 erg/g/s. Only 6.616733e-8 survives, a 99.9799 percent
cancellation. The two face Pi values are -9.42367 and +8.51038. Reconstructing
the code's source over the profile agrees to 6.45e-16 relative to the sum
of the absolute contributions. The averaging mechanism is demonstrably
active in the tail; a controlled replacement is still required to establish
how much of the stellar pattern it causes.

A periodic 32-zone stencil check gives exactly zero alternating response
with the current two averages and a nonzero response when the moments are
colocated. Its one-wavelength smooth response agrees with cos(pi/32)^2 to
roundoff. `output/review/rsp3_modes_profile_20260920/audit_stable_flux.py`
and `stable_flux_coupling.json` preserve the check and profile measurements.

The same profile has weak Pi/Phi turnover decay in these layers. At face
295, the local Phi decay timescale is 4727 days, compared with a local
inverse buoyancy frequency of 2911 seconds. Radiative decay is disabled by
RSP2_alfar=0. The Phi term can therefore continue driving Pi. Enhancing only
energy dissipation does not correct the averaging and, with an unchanged
moment length, reducing w also reduces the existing Pi/Phi decay rates.
The published common-length prescription affects more than energy alone.

This strengthens the recommendation for colocated face energy, Pi and Phi,
with conservative cell energy projection and a consistent local closure.
It does not establish that collocation alone removes the entire tail or
that physical time-dependent heat flux must vanish in every stable layer.

### How Phi is obtained and can acquire spatial structure

Phi is an independent face unknown, the full unresolved entropy variance.
It is not the square of a difference between neighboring cells' mean
entropies. `hydro_rsp2:do1_rsp2_moment_eqns` evolves it with

```math
 \frac{\Phi_k-\Phi_{k,\mathrm{start}}}{\Delta t}
 =2\left(-\frac{\partial s}{\partial r}\right)_k\Pi_k
 -\left[4\sqrt{\frac23}\,\texttt{RSP2\_alfa\_phi}
                 \frac{w_{\mathrm{face},k}}{\Lambda_k}
        +\frac{2}{\tau_{{\rm rad},k}}\right]\Phi_k.
```

The active row uses the new coupled state on its right-hand side. Its
entropy driving comes from `hydro_gradient_support:get_rsp2_thermal_gradient`.
The normal temperature-row paths express it through Y and the matching
pressure, dynamical and composition corrections. The `constant_L` path,
which lacks that temperature row, uses adjacent T and P differences.
Spatial differences enter this driving and other hydro coefficients, but
there is no finite-difference estimate of the variance itself. Its profile
column reads the stored Phi directly.

Spatially alternating Pi can therefore drive spatially uneven production
or destruction of Phi, even with smooth entropy driving. At alfar=0 and
small reconstructed w, its local damping is weak; no nonlocal Phi transport
currently redistributes that structure. The demonstrated cell/face loss of
kinetic-energy response can participate in this coupled pattern. These are
mechanisms supported by the equations and profile audit, not a proof of
the initial seed or the complete cause of every Phi feature.

On initial activation `init_rsp2_moments` seeds active convection with
Phi=1.5*Pi^2/e_face after preserving Lc; quiet stable faces start at zero.
That ratio can inherit sharp input structure and is not an equilibrium
solve. Restart restores Phi, and remesh interpolates it with the existing
monotonicity/positivity guards. Neither path computes it from squared mean
entropy differences. A separate local-closure issue remains: the source
can request negative Phi with nonzero Pi at Phi=0, as documented in the
covariance audit. Changing the spatial gradient estimator cannot alone
resolve that demonstrated homogeneous failure.

## Layout comparison

| Layout | Local moment consistency | Thermal coupling | Future transport and remesh | RSP2 energy accounting |
| --- | --- | --- | --- | --- |
| Cell e_t; face Pi, Phi, as now | Face energy production differs from the production in the local covariance proof | Direct face Pi | Requires compatible coupling between two grids | Retained |
| Cell e_t, Phi; face Pi | Still separates the covariance from its diagonal moments; averaging is insufficient | Direct face Pi | Still requires a proof across two grids | Retained |
| Cell e_t, Pi, Phi with simple averaging | Local consistency can be established | Extra alternating thermal mode in the constant-coefficient test | Common cell transport and overlap remap | Retained |
| Face e_t, Pi, Phi | Local consistency can be established | Direct face Pi and compact thermal response | Common transport and overlap remap on face control volumes | Must be rederived |

Face control volumes extend between adjacent gas cell centers, with half
volumes at the physical boundaries. Face storage does not prohibit
conservative finite-volume transport or require splines. Point values,
volume averages and quadrature weights must be defined consistently on
unequal meshes. Storage location alone does not establish accuracy.

## Heat equation check

Freeze the local thermodynamic coefficients and turbulent energy. Examine
the linear heat and entropy-flux subsystem, with q=2e_t/3 and positive
flux decay rate gamma:

```math
\partial_t\delta s=-\partial_x\delta\Pi,
\qquad
\partial_t\delta\Pi=-q\,\partial_x\delta s-\gamma\delta\Pi.
```

Constant geometric and thermodynamic factors are absorbed in the units.
This is not the full stellar linearization. In the slow limit the heat
diffusivity is q/gamma. With Pi on faces, adjacent entropy values give

```math
\partial_t\delta s_i=\frac{q}{\gamma}
 \frac{\delta s_{i+1}-2\delta s_i+\delta s_{i-1}}{\Delta x^2}.
```

Averaging face gradients into cells and cell Pi back to faces gives instead

```math
\partial_t\delta s_i=\frac{q}{\gamma}
 \frac{\delta s_{i+2}-2\delta s_i+\delta s_{i-2}}{4\Delta x^2}.
```

The latter leaves delta_s_i=(-1)^i unchanged. On a periodic 32-cell grid,
face Pi gives only the physical constant null mode; the simple cell scheme
has two. The extra mode is also present in the time-dependent system.

With q/gamma=1 and radiative diffusivity 0.001, the alternating mode has
decay rate 4.004 with face Pi and 0.004 with averaged cell Pi at unit spacing.
Radiation removes exact singularity but does not restore the missing
convective response. This is a controlled test, not the Cepheid Jacobian.

The time-dependent derivative has Fourier magnitude 2*sin(theta/2)/dx
for face Pi and sin(theta)/dx for the simple cell scheme. At one wavelength
across 32 cells, these are 0.998394 and 0.993587 times the continuum value.
Both are second order for smooth waves. Face Pi has the more accurate wave
response for these stencils. This does not predict a stellar growth rate.

Pulsations require accurate phase, growth and work integrals. Evolution
requires the adjacent-cell thermal response when timesteps exceed turnover
times. Both matter without changing the physical luminosity time weights.
The compact diffusion-limit principle has precedent in
[Kupper, Frank and Jin](https://arxiv.org/abs/1501.02180); their linear
transport method is not a validation of RSP3.

## Nonlocal moments and time discretization

For our full Phi variance and isotropic velocity closure, admissibility is

```math
e_t\geq0,\qquad\Phi\geq0,\qquad\Pi^2\leq(2/3)e_t\Phi.
```

The matrix with entries (2e_t/3, Pi; Pi, Phi) must be positive semidefinite.
A common positive finite-volume diffusion operator mixes these matrices
with nonnegative weights, on either a cell or a face grid. A common positive
overlap remap has the same property and conserves the integrated moments.
Transport/remap e_t, then recover w; do not average w to conserve energy.

Energy-only transport, or arbitrary independent diffusion coefficients,
does not give this guarantee. The earlier energy-export counterexample is
independent of location. Common diffusion is a sufficient candidate, not
a calibration of the physical nonlocal terms. Adding Pi/Phi transport
changes the earlier local-only scope and needs an explicit decision.

Nonlocal transport of Pi is distinct from Lc=4*pi*r^2*rho*T*Pi. The latter
transports mean enthalpy; the former transports the velocity/entropy
covariance. For face moments, their nonlocal transport fluxes lie between
moment locations, at gas cell centers. They need not coincide with Lc.

The stellar literature does not prove that nonlocal terms cure a bad
local closure. [Braun et al.](https://arxiv.org/html/2604.06151v1) show that
local and nonlocal closures can substantially change mean stratification.
Admissibility, calibration and realistic temperature profiles are separate
requirements.

A common spatial operator alone is insufficient. In a three-cell example
at dt=10, Crank-Nicolson maps variance (0,1,0) to central variance -0.25.
Backward Euler has positive mixing weights for the same diffusion operator.
This is not a proposal to make total luminosity fully implicit. RSP2
currently time centers Lt while Pi/Phi source rows use the new state.
Future moment transport needs a consistent discrete time update and
separate pulsation accuracy checks. Implicitness alone proves neither
admissibility nor accuracy.

## Conservative face energy

An all-face formulation must not substitute face w into existing cell
expressions. Conservative storage and transport can be constructed.
For this algebra, indices increase in the transport direction; MESA's
inward indices and outward luminosity require the corresponding signs.
Cell i lies between faces i and i+1. Define

```math
e_{t,\mathrm{cell},i}=\tfrac12(e_{t,i}+e_{t,i+1}),\qquad
\Delta m_{\mathrm{face},i}=\tfrac12(\Delta m_{i-1}+\Delta m_i).
```

The end faces have half-cell masses. Exactly,

```math
\sum_i\Delta m_i e_{t,\mathrm{cell},i}
 =\sum_j\Delta m_{\mathrm{face},j}e_{t,j}.
```

If face energy evolves by conservative fluxes Lt_dual at cell centers,
its cell counterpart uses the face flux

```math
L_{t,j}=
 \frac{\Delta m_j L^{\rm dual}_{t,j-1/2}
       +\Delta m_{j-1}L^{\rm dual}_{t,j+1/2}}
      {\Delta m_{j-1}+\Delta m_j}.
```

Its divergence gives the same cell energy derivative as averaging the
two face derivatives. The unequal-mass standalone check agrees within
4.5e-16. This proves storage and diffusion identities at fixed cell masses,
not the entire hydro energy identity.

Pressure work, Eq/Uq, boundary suppression, u/v forms, both thermal energy
equations and time centering still need compatible derivations. A surface
Dirichlet condition cannot silently discard a half-volume's energy.
Remeshing must preserve thermal plus mechanical energy on the new volumes.

TDC already supplies a face-velocity precedent in `hydro_energy.f90`,
`setup_d_turbulent_energy_dt`. Its normalization is different: its current
cell energy is 0.75*(vc_face^2+vc_next_face^2). This is not permission to
merge TDC/RSP2 hydro work routines or assume their identities transfer.

## Source map and verification required

Relevant current paths:

- `star/private/hydro_rsp2.f90`: `rsp2_moment_rhs`, `rsp2_moment_source`,
  `do1_turbulent_energy_eqn`, `compute_Lc_terms`, `compute_Lt`,
  `rsp2_flux_residual` and boundary rules.
- `star/private/hydro_gradient_support.f90`: `get_rsp2_thermal_gradient`
  maps each temperature row and dynamical/composition terms to face entropy
  driving. Preserve that mapping and its AD dependencies.
- `star/private/hydro_energy.f90`: cell turbulent storage.
- `star/private/star_utils.f90`: energy diagnostics and `calc_Ptrb_ad_tw`
  currently interpret RSP2 w as a cell quantity.
- `star/private/mesh_adjust.f90`: `do_etrb` conserves cell dm*w^2;
  Pi/Phi use separate face interpolation.
- `star/private/adjust_mesh_split_merge.f90` and
  `star/private/tdc_hydro_support.f90`: split/merge and envelope remesh.
- Allocation, `set_flags`, AD wrappers, photo/model I/O, profiles, mixing
  and `star_LNA*`: location-dependent consumers.

Before enabling the proposed layout:

1. Finish the local closure and discrete admissibility derivation, including
   zero moments, compression and Newton trials.
2. Derive work, viscous heating and transport using the same cell energy
   projection. Preserve total L time centering unless separately justified.
3. Define boundary contributions and common remap in all three mesh paths.
   Check strongly unequal masses and actual boundary conditions.
4. Specify restart conversion. Old cell w and new face w cannot share a
   photo interpretation. The failing photo already contains incompatible
   moments; relocation alone cannot repair them.
5. Preserve one-equation RSP2 behavior and linearize the selected RSP3
   equations consistently in LNA. Update the manuscript only for the
   implementation actually selected.
6. Check AD partials, steady thermal response and resolved wave dispersion.
   Then compare periods, growth, energy error, retries and timesteps over
   several stellar pulsation cycles.
7. Test stellar evolution with timesteps longer than turnover times,
   changing convection boundaries and composition gradients. Test both
   local and nonlocal limits without using radiation to hide missing
   convective coupling.

The result is a better-supported design recommendation. Neither the current
closure nor an unimplemented replacement is declared robust by this audit.

## Reproducible checks

`output/review/rsp3_restart_115000_20260920/check_layouts.py` writes
`layout_checks.json`. It checks the heat stencils and their time-dependent
systems, common diffusion, a time-centering counterexample, unequal-mass
energy projection and overlap remapping. All assertions pass. Common
diffusion integral errors are below 5e-15 relative; remap covariance is
nonnegative to roundoff. These are mathematical checks, not stellar runs.
The earlier two stellar restarts remain evidence against the current
closure only.
