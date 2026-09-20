# RSP3 moment and enthalpy closure audit

Date: 2026-09-20. Source comparison only. No change to the physical equations,
controls or test case. No MESA compilation or execution.

Sources checked:

- Kovacs, Szabo and Nuspl (2026), [arXiv v1](https://arxiv.org/html/2601.04931v1),
  equations (3), (5)--(15), sections 3.3--3.5 and appendices A and C.
  Local PDF: `references/kovacs_szabo_nuspl_2026_three_equation_tdc.pdf`.
  Printed pages 3, 7 and 8 were inspected. The publisher version and an
  author correction have not been verified.
- Braun, Ahlborn, Kupka and Weiss (2026),
  [section 2](https://arxiv.org/html/2604.06151v1#S2), including the local
  and nonlocal alternatives for the third order moments.

## Definitions and implemented equations

The code uses `Pi = <v_r' s'>`, `Phi = <s'^2>` and cell energy `w**2`.
Kovacs and Braun use half the entropy variance. In this note only,
`Phi_K = Phi/2` distinguishes their normalization. It is not a new
code variable.

In `star/private/hydro_rsp2.f90`, `rsp2_moment_rhs` implements the
following face equations. The displayed derivatives stand for the
existing discrete face gradient and velocity strain. The face value
`w_face**2` is reconstructed from adjacent cell energies.

```math
 \frac{D\Pi}{Dt}=
 -\frac23 w_{\rm face}^2\frac{\partial s}{\partial r}
 -\frac{\delta\Phi}{\rho c_p}\frac{\partial P}{\partial r}
 -\left[
 6\sqrt{2/3}\,{\tt RSP2\_alfa\_pi}\frac{w_{\rm face}}{\Lambda}
 +\frac1{\tau_{\rm rad}}+\frac{\partial v}{\partial r}
 \right]\Pi,
```

```math
 \frac{D\Phi}{Dt}=-2\Pi\frac{\partial s}{\partial r}
 -\left[
 4\sqrt{2/3}\,{\tt RSP2\_alfa\_phi}\frac{w_{\rm face}}{\Lambda}
 +\frac2{\tau_{\rm rad}}
 \right]\Phi.
```

Here `delta = ChiT/ChiRho`, and

```math
 \tau_{\rm rad}^{-1}=
 \frac{48\,{\tt RSP2\_alfar}^{,2}\sigma T^3}
      {c_p\kappa\rho^2\Lambda^2}.
```

`x_ALFAPI` and `x_ALFAPHI` hold the local calibration constants.
The code variable `entropy_gradient` supplies the negative entropy
gradient represented above, with the existing thermal and composition
mapping. `rsp2_buoyancy_face` supplies the pressure gradient coefficient.
`compute_Lc_terms` implements

```math
 L_c=4\pi r^2\rho_{\rm face}T_{\rm face}\Pi_{\rm face}.
```

These formulas establish that the implementation already has local
moment decay and radiative cooling. They do not establish that its
physical closure is complete or realizable.

## Transport and additional dissipation

Braun's local closure replaces a transport divergence by
`alpha_a*w*a/Lambda`. This is the origin of our two local decay terms.
A rate with units `a/time` is not a spatial flux. Our local equations
do not supply the flux needed in a higher order enthalpy expression.
They also do not prove that this flux is physically zero. Recovering
one from the local rate would require an additional spatial construction
and boundary conditions.

Kovacs equation (7) instead uses the actual variance flux

```math
 \mathcal F_{\Phi_K}=\frac\rho2\langle v_r's'^2\rangle,
 \qquad
 F_c=\rho T\Pi+\alpha_c\frac{c_{pT}T}{c_p}\mathcal F_{\Phi_K}.
```

The paper treats this as an ionization correction. Its epsilon terms
include ordinary radiation and extra viscous decay; the separate
tau_kappa terms describe opacity fluctuations. Expanding those terms
and converting to our full variance gives

```math
 \left.\dot\Pi\right|_{\rm loss}=
 -\left[\frac{4-\kappa_T}{\tau_r}
       +\alpha_{d,\Pi}\frac{w}{\Lambda}\right]\Pi,
 \qquad
 \left.\dot\Phi\right|_{\rm loss}=
 -\left[\frac{2(4-\kappa_T)}{\tau_r}
       +2\alpha_{d,\Phi}\frac{w}{\Lambda}\right]\Phi.
```

These alpha_d coefficients are not our alpha controls. Relative to our
source, the proposed extra viscous terms would be additions to the
existing local transport closure. The opacity terms would require a
consistent mean energy treatment as well. Neither extension is present.
With `RSP2_alfar=0`, our inverse cooling time is zero. Merely multiplying
it by an opacity derivative would leave it zero.

## Printed inconsistencies

The half variance definition requires twice the printed Phi buoyancy
term in Kovacs equation (10). Equation (13) and appendix A also differ
by two in viscous variance decay. Appendix A's opacity driving inequality
disagrees with its displayed timescale. These are internal discrepancies,
not an author confirmed erratum.

An independent normalization check is direct. At constant pressure,
`rho'/rho = -delta*s'/cp`. Multiplying the fluctuating buoyancy
acceleration by `s'` gives

```math
 -\frac{\delta\langle s'^2\rangle}{\rho c_p}\partial_rP
 =-\frac{2\delta\Phi_K}{\rho c_p}\partial_rP
 =-\frac{\delta\Phi}{\rho c_p}\partial_rP.
```

Our source has the last form. Halving it to reproduce the printed
equation would introduce a normalization error. This check is distinct
from the proposed pressure redistribution and covariance correction in
`rsp3_covariance_closure.md`.

## Independent enthalpy check before adopting equation (7)

For a thermodynamically consistent equilibrium EOS at fixed pressure
and nuclear composition, `dh = T ds`. Consequently

```math
 h_s=T,\qquad h_{ss}=T/c_p,
 \qquad
 h'=Ts'+\frac{T}{2c_p}(s'^2-\langle s'^2\rangle)+\cdots .
```

Keeping the same leading density approximation as our current flux gives

```math
 \rho\langle v_r'h'\rangle
 =\rho T\Pi+\frac{T}{c_p}\mathcal F_{\Phi_K}+\cdots .
```

This coefficient has no `cpT` factor. Expanding enthalpy in temperature
and converting temperature to entropy consistently through second order
gives the same result:

```math
 T'=\frac{T}{c_p}s'
 +\frac{T(1-c_{pT})}{2c_p^2}(s'^2-\langle s'^2\rangle)+\cdots .
```

The contribution from this quadratic conversion cancels the `cpT`
dependence from `h_TT`. This is a restricted thermodynamic check, not
a complete compressible mean flux derivation. It does show why copying
the ionization coefficient into MESA requires resolving the averaging
and EOS assumptions first. Independent ionization or composition
fluctuations would require a separate closure.

## Consequence for the present diagnosis

The source audit establishes which extensions are absent. It does not
establish that any of them caused the current broad turbulent tail.
The saved profile energy budget in `rsp2_three_equation_implementation.md`
still identifies continued shear production and weak local braking.
Changing the enthalpy flux, adding moment dissipation and changing stable
layer energy dissipation are distinct physical changes. None is applied
by this audit.

Before implementing the Felix Ahlborn dissipation prescription, derive
its mapping to our energy and mixing length conventions, its zero-w
limit, gas energy exchange, and LNA derivatives. Keep that decision
separate from the pending covariance and face layout work.

## Local relaxation versus viscous decay

Question to revisit: Are our local Pi/Phi relaxation terms mimicking
Kovacs's viscous terms, and should we replace them?

Short answer: they produce the same mathematical damping response in
the local equations, but approximate different physical processes.
Our terms stand in for transport between layers. His additional terms
represent viscous destruction of correlations. Our local replacement
does not actually transfer anything to a neighboring layer. Calling the
terms physically distinct must not obscure their identical local
mathematical form.

Both prescriptions have the same algebraic dependence on `w*Pi/Lambda`
and `w*Phi/Lambda`. In the present local equations their separate
coefficients would therefore add. They are distinguished by the terms
they approximate, not by a different response in the local solver.

An actual transport divergence can import or export a moment; its mass
integral reduces to boundary fluxes. The local replacement instead
relaxes the moment toward zero and does not pass it to an adjacent zone.
The explicit viscous model describes local destruction of correlations.
The local transport replacement is consequently a substantive closure
approximation, not a conservative implementation of nonlocal transport.

If both prescriptions were included literally, their turnover terms in
our full variance convention would be

```math
 \left.\dot\Pi\right|_{\rm turnover}=
 -\left[6\sqrt{2/3}\,{\tt RSP2\_alfa\_pi}+\alpha_{d,\Pi}\right]
 \frac{w}{\Lambda}\Pi,
```

```math
 \left.\dot\Phi\right|_{\rm turnover}=
 -\left[4\sqrt{2/3}\,{\tt RSP2\_alfa\_phi}+2\alpha_{d,\Phi}\right]
 \frac{w}{\Lambda}\Phi.
```

| Coefficient of turnover decay | Current unit control defaults | Additional Kovacs term |
| --- | ---: | ---: |
| Pi | 4.898979486 | 0.0217 |
| Phi | 3.265986324 | 0.0434 |

Using the quoted Kovacs coefficients `alpha_d,Pi = alpha_d,Phi = 0.0217`
and his main equation (13), the extra turnover decay coefficients in
our full variance convention would be 0.0217 and 0.0434. Our unit control
defaults give 4.898979486 and 3.265986324. Literal addition increases those
coefficients by about 0.44 and 1.33 percent. Replacement would instead
reduce them by factors of about 226 and 75. These ratios concern only
turnover decay, not total damping including radiation or mean strain.
The printed variance factor discrepancy remains unresolved.

There is no basis to replace the present terms as a bug fix. Choosing
zero spatial transport with only the small viscous terms defines another
local model and abandons the current local calibration. If actual Pi/Phi
transport is implemented later, it replaces the corresponding local
transport approximation; a separately justified viscous term can remain.

Status at this discussion: retain the current equations. Neither the
replacement nor the addition has been implemented. The current terms
are in `rsp2_moment_rhs`, with constants `x_ALFAPI` and `x_ALFAPHI`, in
`star/private/hydro_rsp2.f90`. Recheck that source before describing a
later implementation. This comparison does not establish a remedy for
the instability or excessive mixing, and is separate from the proposed
Ahlborn stable layer energy dissipation change.
