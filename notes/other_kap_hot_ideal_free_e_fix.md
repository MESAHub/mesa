# Opacity-Only Replacement of the Ideal-EOS Electron Component

Date: 2026-08-22

Status: the hot ideal-electron correction is implemented in the case-local
`my_other_kap_get`. Static Fortran lint passes. The low-temperature opacity
continuation documented at the end of this note is also implemented in the
same callback. No compilation or MESA run has been performed for the
continuation.

Scope: the case-local `my_other_kap_get` in both PPISN development cases and
`star/test_suite/ppisn/src/run_star_extras.f90`. This contains the opacity
failure without changing the EOS thermodynamic state. The separate HELM
selector/evaluator mismatch is documented in `notes/helm_bug_fix.md`.

## What is wrong

The relevant end of the EOS hierarchy is

```text
FreeEOS -> OPAL/SCVH -> HELM -> ideal
```

The ideal EOS is the final coverage fallback. It does not calculate an
ionization state. In `eos/private/ideal.f90:ideal_eos`, MESA sets

```fortran
! No electrons, so extreme negative chemical potential
etaele = -1d99
xnefer = 1d-20
```

Here `xnefer` is an electron number density, not the dimensionless electrons
per baryon passed to opacity. `eos/private/skye_thermodynamics.f90:pack_for_export`
converts it to

\[
\ell_I
=\ln f_I
=\ln\left(\frac{10^{-20}}{N_A\rho}\right).
\]

Therefore the raw ideal-component derivatives at fixed composition are

\[
\frac{\partial\ell_I}{\partial\ln\rho}=-1,
\qquad
\frac{\partial\ell_I}{\partial\ln T}=0.
\]

The resulting near-zero `free_e` is not evidence for physical recombination.
It is the neutral value attached to a deliberately minimal coverage fallback.

`eos/private/eosdt_support.f90:Do_Blend` blends every EOS result, including
`lnfree_e`, with the same quintic-smoothed EOS weights. Passing the final
blended `lnfree_e` into `kap/private/kap_eval.f90:Compton_Opacity` then gives

\[
\kappa_{\rm C}\propto\exp({\tt lnfree\_e}),
\]

so even a modest ideal fraction can collapse the Compton opacity by many
orders of magnitude.

The correction should not add a separate density blend or an ionization
model. It should replace only the ideal EOS contribution to `lnfree_e`, using
the same ideal fraction and the same FreeEOS/HELM-to-ideal blend already
calculated by the EOS.

## The EOS boundaries are diagnostic, not new callback conditions

The diagonal FreeEOS boundary uses

\[
\log Q=\log_{10}\rho-2\log_{10}T+12.
\]

In `eos/private/eosdt_eval.f90:Get_FreeEOS_alfa`, the default low-`Q`
transition is

\[
\begin{aligned}
\log Q &< -10.0: && f_{\rm FreeEOS}=0,\\
-10.0\leq\log Q&<-9.9:
   && f_{\rm FreeEOS}=\frac{\log Q+10.0}{0.1},\\
\log Q&\geq-9.9: && f_{\rm FreeEOS}=1,
\end{aligned}
\]

when the other FreeEOS availability factors are unity. The corresponding
lines are

\[
\log_{10}\rho=2\log_{10}T-22.0
\]

and

\[
\log_{10}\rho=2\log_{10}T-21.9.
\]

The actual lower HELM table coordinate is instead

\[
\log_{10}(\rho Y_e)=-12,
\qquad
Y_e=\sum_i X_i\frac{Z_i}{A_i}.
\]

The EOS-region plot for `X = 0`, `Z = 1` shows that the coverage boundaries
extend to substantially lower temperatures than the model-11550 failure.
A high cutoff near `logT = 4.5` would therefore miss part of the coverage
fallback.

Neither `logQ` nor `rho*Ye` should be reconstructed in `other_kap`. The final
EOS result already stores

```fortran
a = s% eos_frac_ideal(k)
```

after all selectors, nested blends, composition cuts, and converted EOS
evaluation failures have been applied. This is both the coverage indicator
and the correct weight for replacing the ideal electron component.

## Replacement ideal component

Calculate the composition electron abundance

\[
Y_e=\sum_i X_i\frac{Z_i}{A_i}=\frac{\bar Z}{\bar A}
\]

with `chem_lib:basic_composition_info` from the callback's `species`,
`chem_id`, and `xa` arguments.

The scope floor remains `log10_T = 3.7`. Above that floor, use the
opacity-only ideal-component electron abundance

\[
f_{I,{\rm kap}}=Y_e,
\]

where

\[
f_I=\frac{10^{-20}}{N_A\rho}.
\]

Then

\[
\ell_{I,{\rm kap}}=\ln f_{I,{\rm kap}}
\]

replaces only the ideal component's original

\[
\ell_I=\ln f_I.
\]

Define

\[
\Delta_I=\ell_{I,{\rm kap}}-\ell_I.
\]

This is not a new ionization calculation. It is the opacity convention used
for the otherwise undefined electron content of the ideal coverage fallback.

## Full expression using the existing EOS blend

Let

\[
\ell_0={\tt lnfree\_e}
\]

be the final blended EOS value received by `other_kap`, and let

\[
a=f_{\rm ideal}={\tt s\%eos\_frac\_ideal(k)}.
\]

Because `Do_Blend` combines `lnfree_e` linearly with the EOS weights, replacing
only the ideal component gives the exact corrected blend

\[
\boxed{
\ell_{\rm kap}=\ell_0+a\Delta_I
}
\]

or, equivalently,

\[
\boxed{
f_{e,{\rm kap}}
=f_0\exp(a\Delta_I),
\qquad
f_0=\exp(\ell_0)
}.
\]

This has the desired limits:

\[
\begin{array}{lll}
a=0: & \ell_{\rm kap}=\ell_0,
   & \text{no ideal contribution},\\
0<a<1: & \ell_{\rm kap}=\ell_0+a\Delta_I,
   & \text{same EOS transition width},\\
a=1: & \ell_{\rm kap}=\ell_{I,{\rm kap}},
   & \text{pure ideal coverage fallback}.
\end{array}
\]

There is no `max(f_0,Y_e)` and no separate density smoothstep. The non-ideal
FreeEOS or HELM component remains exactly as calculated; only the weighted
ideal component is exchanged.

## Exact logarithmic derivatives

Let

\[
a_\rho=\frac{\partial a}{\partial\ln\rho},
\qquad
a_T=\frac{\partial a}{\partial\ln T},
\]

and denote the input derivatives by

\[
\ell_{0,\rho}
=\frac{\partial\ell_0}{\partial\ln\rho},
\qquad
\ell_{0,T}
=\frac{\partial\ell_0}{\partial\ln T}.
\]

At fixed composition, \(\ell_{I,{\rm kap}}=\ln Y_e\) has zero density and
temperature derivatives. Since \(\ell_{I,\rho}=-1\) and
\(\ell_{I,T}=0\),

\[
\boxed{
\Delta_{I,\rho}=1
}
\]

and

\[
\boxed{
\Delta_{I,T}=0
}.
\]

The exact derivatives corresponding to the corrected value are therefore

\[
\boxed{
\frac{\partial\ell_{\rm kap}}{\partial\ln\rho}
=\ell_{0,\rho}+a_\rho\Delta_I+a
}
\]

and

\[
\boxed{
\frac{\partial\ell_{\rm kap}}{\partial\ln T}
=\ell_{0,T}+a_T\Delta_I
}.
\]

`a`, `a_rho`, and `a_T` supply the existing FreeEOS/HELM-to-ideal transition
geometry. The `a_rho*Delta_I` and `a_T*Delta_I` terms are essential because
they describe movement through that EOS blend.

## Required access to the ideal-fraction derivatives

The EOS already calculates and stores the required derivatives. `Do_Blend`
returns them as

```fortran
d_dlnd(i_frac_ideal)
d_dlnT(i_frac_ideal)
```

using the same quintic-smoothed blend as `a`. `do_eos_for_cell` passes the
complete derivative vectors directly into the existing step-work arrays:

```fortran
s% d_eos_dlnd(:,k)
s% d_eos_dlnT(:,k)
```

Therefore the callback can use

```fortran
s% d_eos_dlnd(i_frac_ideal,k)
s% d_eos_dlnT(i_frac_ideal,k)
```

directly. No new `star_info` fields, allocation, copying, remeshing, or photo
state are required. `store_eos_for_cell` separately stores the scalar
`s% eos_frac_ideal(k)`, but no additional scalar derivative copies are needed.

## Callback structure

Use local variables for

```fortran
lnfree_e_for_kap
d_lnfree_e_for_kap_dlnRho
d_lnfree_e_for_kap_dlnT
```

initialized from the callback inputs. For a valid cell with
`s% eos_frac_ideal(k) > 0` and `log10_T > 3.7`, calculate `Y_e`, the raw ideal
component, its opacity replacement, and the corrected blended value and
derivatives above. Then make the normal call:

```fortran
call kap_get( &
   handle, species, chem_id, net_iso, xa, &
   log10_rho, log10_T, &
   lnfree_e_for_kap, &
   d_lnfree_e_for_kap_dlnRho, d_lnfree_e_for_kap_dlnT, &
   eta, d_eta_dlnRho, d_eta_dlnT, &
   kap_fracs, kap, dln_kap_dlnRho, dln_kap_dlnT, &
   dln_kap_dxa, ierr)
```

Use the callback's `handle` rather than `s%kap_handle`. Do not overwrite
`s%lnfree_e`, `s%eta`, or any EOS thermodynamic array. The existing unbound
`k == 1` electron-scattering override retains precedence.

The correction has no saved state and depends only on the current trial cell.
It therefore behaves the same on fresh starts, photo restarts, retries, solver
iterations, remeshing, and OpenMP cell evaluation.

The implementation is in the `my_other_kap_get` callbacks under
`star/dev_cases_test_TDC/ppisn_models/dev_TDC_ppisn_from_binary`,
`star/dev_cases_test_TDC/ppisn_models/dev_TDC_single_star_ppisn`, and
`star/test_suite/ppisn`. It also initializes `kap_fracs` before the existing
surface-opacity early return and passes the callback's `handle` into
`kap_get`.

## Model-11550 evidence

The two coverage crossings are visible in profile 141:

| Zone | `logT` | `logQ` | `logRhoYe` | low-`Q` FreeEOS factor | `free_e` |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 30 | 4.48904 | -9.66033 | -12.98287 | 1.0000 | 0.4988 |
| 31 | 4.68935 | -10.55286 | -13.47478 | 0.0000 | \(2.48\times10^{-31}\) |
| 132 | 5.09909 | -9.99435 | -12.09682 | 0.0565 | \(1.17\times10^{-32}\) |
| 137 | 5.12225 | -9.96366 | -12.01982 | 0.3634 | \(1.20\times10^{-24}\) |
| 138 | 5.17195 | -10.04325 | -12.00000 | 0.0000 | \(8.31\times10^{-33}\) |
| 139 | 5.20405 | -10.09921 | -11.99176 | 0.0000 | 0.50044 |

Zone 31 crosses the low-`Q` FreeEOS edge while HELM is unavailable. Zone 139
crosses back into HELM. The intermediate ideal fractions already define the
finite-width density transition required by the opacity correction.

For zone 138, `logT > 3.7`. If the cell is pure ideal fallback, the opacity
input becomes approximately `Y_e ~= 0.5`. In the
nonrelativistic limit this restores

\[
\kappa_{\rm es}
\simeq\frac{\sigma_T}{m_u}Y_e
\simeq0.20\ {\rm cm^2\,g^{-1}},
\]

consistent with the neighboring HELM result.

## Validation plan

The user compiles and runs MESA. After implementation, inspect:

1. model 11550 zones 25 through 145, especially zones 30/31 and 138/139;
2. `eos_frac_ideal` and its two stored derivatives through both coverage
   crossings;
3. `log_opacity`, `kap_frac_Compton`, `dlnkap_dlnRho`, and `dlnkap_dlnT`
   across the same crossings;
4. cells below `logT = 3.7`, which must remain unchanged;
5. cells with `eos_frac_ideal = 0`, which must be exactly unchanged;
6. the existing unbound surface-cell opacity override and all `kap_fracs`
   outputs;
7. a photo restart, retry, and remesh through the affected layer;
8. `solver_test_kap_partials` below and above the EOS blend.

## Low-temperature opacity-table continuation

Status: implemented in the case-local `my_other_kap_get`; static Fortran lint
passes, but the change has not been compiled or run.

This is a separate numerical safeguard from the hot ideal-electron correction
above. The ideal-electron correction repairs the `free_e` input to opacity for
`log10_T > 3.7` when the EOS uses its ideal coverage fallback. The continuation
below instead prevents `kap_get` from being evaluated outside its low-
temperature table coverage. It does not use `eos_frac_ideal` and does not alter
the EOS state.

The late failure at model 17643 repeatedly places zone 163 at approximately

\[
\log_{10}\rho=-9.7735,
\qquad
\log_{10}T=2.7000,
\]

which is the opacity table's lower temperature edge. Retrying to progressively
smaller timesteps does not move the trial state away from that boundary. The
purpose of this continuation is to give the solver a finite opacity and a
continuous first temperature derivative when a trial cell crosses the table
edge.

### Definition

Let

\[
x=\log_{10}T,
\qquad
x_0=2.7001,
\qquad
x_1=x_0+\Delta x,
\qquad
\Delta x=0.1.
\]

Here `x0` is slightly inside the valid table rather than exactly on its
boundary. At fixed density, composition, and opacity auxiliary inputs, define

\[
\ell(x)=\ln\kappa(\rho,T,X_i),
\qquad
\ell_0=\ln\kappa(\rho,T_0,X_i),
\qquad
T_0=10^{x_0}\ {\rm K}.
\]

Within the transition interval, use

\[
u=\frac{x-x_0}{\Delta x},
\qquad
w(u)=3u^2-2u^3.
\]

The extended log opacity is

\[
\boxed{
\ell_{\rm ext}(x)=
\begin{cases}
\ell_0, & x\leq x_0,\\[3pt]
(1-w)\ell_0+w\ell(x), & x_0<x<x_1,\\[3pt]
\ell(x), & x\geq x_1.
\end{cases}
}
\]

Blending `ln(kap)` rather than `kap` guarantees a positive opacity and avoids
allowing the larger of two opacity values to dominate merely because the blend
is performed in linear units.

### What it does as temperature changes

Below `x0`, the opacity is flat as a function of temperature at fixed density
and composition:

\[
\kappa_{\rm ext}(\rho,T,X_i)
=\kappa(\rho,T_0,X_i),
\qquad x\leq x_0.
\]

It therefore neither decreases nor increases with further cooling below the
table edge. This is not one global constant opacity: the boundary value is
re-evaluated for the cell's current density and composition, so different cells
can have different continued opacities.

Between `x0 = 2.7001` and `x1 = 2.8001`, the opacity joins smoothly onto the
normal table. It does not necessarily "smoothly decrease." Its direction is
set by the actual opacity table over that interval: it rises with temperature
if the table value rises, and falls if the table value falls. Above `x1`, the
normal opacity is exactly unchanged.

Along an evolving stellar profile, density and composition generally change at
the same time as temperature. The continued opacity can therefore change from
one zone or timestep to the next even while that zone remains below `x0`; only
its partial derivative with respect to temperature at fixed density and
composition is zero there.

### Consistent logarithmic derivatives

Because `w` depends only on temperature, the density derivative in the blend
interval is

\[
\boxed{
\frac{\partial\ell_{\rm ext}}{\partial\ln\rho}
=(1-w)\frac{\partial\ell_0}{\partial\ln\rho}
+w\frac{\partial\ell}{\partial\ln\rho}
}.
\]

Since

\[
\frac{dw}{d\ln T}
=\frac{6u(1-u)}{\ln(10)\,\Delta x},
\]

the temperature derivative is

\[
\boxed{
\frac{\partial\ell_{\rm ext}}{\partial\ln T}
=w\frac{\partial\ell}{\partial\ln T}
+\frac{6u(1-u)}{\ln(10)\,\Delta x}
  \left(\ell-\ell_0\right)
}.
\]

The composition derivatives are blended in the same way as the density
derivative:

\[
\boxed{
\frac{\partial\ell_{\rm ext}}{\partial X_i}
=(1-w)\frac{\partial\ell_0}{\partial X_i}
+w\frac{\partial\ell}{\partial X_i}
}.
\]

Below `x0`, set

\[
\frac{\partial\ell_{\rm ext}}{\partial\ln T}=0
\]

while retaining the density and composition derivatives returned by `kap_get`
at `x0`. At both endpoints, `dw/dlnT = 0`. In addition,
`ell(x0) = ell0`, so the opacity and its first temperature derivative are
continuous at `x0`; at `x1`, the extended value and derivative equal the normal
table value and derivative. The construction is therefore `C1`, which is the
continuity needed by the Newton solve. A quintic smoothstep could also make the
second derivative continuous, but is not required for the opacity partials
used by the solver.

`kap_fracs` are diagnostics rather than logarithmic quantities. Below `x0`,
use the fractions from the boundary evaluation; within the transition, blend
the two fraction vectors linearly with the same `w`.

### Placement in `my_other_kap_get`

The continuation is applied only around the normal `kap_get` paths in those
three `my_other_kap_get` callbacks. The existing unbound `k == 1` electron-
scattering override must keep precedence.

The implemented evaluation order is:

1. calculate the already implemented hot ideal-electron replacement inputs;
2. if the surface electron-scattering override applies, return it unchanged;
3. for `x >= x1`, make the existing single `kap_get` call at `x`;
4. for `x0 < x < x1`, call `kap_get` at both `x0` and `x`, then combine the
   values and derivatives using the equations above;
5. for `x <= x0`, call `kap_get` only at `x0`, return its opacity, density and
   composition derivatives, and set `dln_kap_dlnT = 0`.

The boundary call must use `x0`, so `kap_get` is never asked to extrapolate
below its table. If the call at `x0` itself reports an error, propagate that
`ierr`; do not hide an unrelated density, composition, or table-availability
failure.

The opacity callback also receives `lnfree_e`, `eta`, and their derivatives
from the EOS evaluation at the trial cell. Strictly, those auxiliary inputs are
not re-evaluated at `T0`. At these cool temperatures the failing path is the
radiative opacity table rather than the hot Compton contribution for which
those inputs matter. This is therefore deliberately a numerical table-edge
continuation, not a new low-temperature microphysics model.

### Expected effect and limitation

This removes the abrupt `kap_get` domain failure and removes a kink in the
opacity supplied to the solver. It also prevents the opacity from continuing
to collapse solely because a trial cell cools below the table base. It does not
force the opacity upward or downward above the table base, and it does not by
itself prevent the hydrodynamic temperature variable from crossing MESA's own
minimum allowed `logT`. The run must still be checked to determine whether
making the opacity equation well defined is enough to stop the rapid cooling
and timestep collapse.

Validate the value and all returned partials with a small scan across
`log10_T = 2.65--2.85` at representative densities and compositions from the
failing ejecta. In particular, compare finite-difference partials immediately
below, inside, and immediately above both transition endpoints. The user
compiles and runs MESA.
