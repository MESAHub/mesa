# RSP3 covariance closure correction

Date: 2026-09-20. Status: derivation and standalone checks, not implemented.
The two requested stellar restarts are complete. No production source was
edited; all test input edits were confined to isolated copies. The user authorized deriving
a consistent correction after the variance failure was identified.

The subsequent [placement audit](rsp3_layout_audit.md) compares the current
staggering, cell moments and face moments. Its preferred simple long-term
RSP3 design places all three moments on faces, after rederiving the energy
projection. It identifies a thermal coupling defect in the simple cell
averaging candidate below. Neither layout change is implemented. The local
closure derivation remains applicable to moments evaluated together.

## Run evidence

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

This proposal adds Pi and Phi transport, contrary to the earlier requested
local-only approximation. It requires a scope decision before coding.
It also does not guarantee the absence of a negative mean temperature
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
