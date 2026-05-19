---
title: "EOS Composition Partials: Implementation Map v2"
subtitle: "Current row-scoped EOS dxa, Skye speed path, Brunt, and implicit diffusion state"
date: "2026-05-17"
geometry: margin=0.55in
fontsize: 10pt
colorlinks: true
linkcolor: blue
urlcolor: blue
header-includes:
  - \usepackage{amsmath}
  - \usepackage{amssymb}
  - \usepackage{booktabs}
  - \usepackage{array}
---

This is the current implementation map for the development branch.  It
supersedes `notes/eos_composition_partials_implementation_map.md` as an
architecture snapshot.  The older file is still useful as the original plan.

Latest verification status: prior install/plotter validation exists from
before the latest Skye Coulomb speed refactors.  The latest edits in this map
have had static checks only; MESA compile/model runs are user-owned.

# Current Goal

The production goal is to make EOS composition partials available where the
star solver actually consumes them, without replacing physics with limiters or
finite-difference approximations.

The current production star-solver request is row-scoped:

\[
\left\{
  \frac{\partial \ln P_{\rm gas}}{\partial X_j},
  \frac{\partial \ln E}{\partial X_j}
\right\}_{\rho,T},
\]

with wider rows available for diagnostics and future consumers.

# Primary Files

| Area | Main files |
|---|---|
| EOS composition helper | `eos/private/eos_composition_partials.f90` |
| eosDT row/moment plumbing | `eos/private/eosdt_eval.f90`, `eos/public/eos_lib.f90` |
| Skye EOS partials | `eos/private/skye.f90`, `eos/private/skye_ideal.f90`, `eos/private/skye_coulomb.f90`, `eos/private/skye_composition_ad.f90` |
| Skye Coulomb leaves | `eos/private/skye_coulomb_liquid.f90`, `eos/private/skye_coulomb_solid.f90` |
| FreeEOS/OPAL table mapping | `eos/private/eosdt_eval.f90` |
| Timing diagnostics | `eos/private/eos_timing.f90`, `star/job/run_star_support.f90` |
| Star EOS calls | `star/private/eos_support.f90`, `star/private/micro.f90` |
| Implicit Brunt/Dmix | `star/private/implicit_brunt.f90`, `star/private/implicit_Dmix.f90`, `star/private/hydro_vars.f90`, `star/private/hydro_chem_eqns.f90`, `star/private/mix_info.f90` |

# Composition Gauge

MESA evolves constrained mass fractions:

\[
\sum_i X_i = 1, \qquad \sum_i \delta X_i = 0.
\]

Composition derivatives are therefore defined up to an additive constant:

\[
\sum_i (Q_i + C)\delta X_i = \sum_i Q_i\delta X_i.
\]

The star solver still applies its sink projection in
`star/private/star_solver.f90`:

\[
D_j^{\rm sink} Q = Q_j - Q_s.
\]

The EOS side computes a consistent constrained isotope-column gauge from
composition moments, then the star-side sink projection can be applied as
usual.

# Moment Derivatives

The shared helper provides the raw composition moment columns:

\[
\bar A = \frac{\sum_i X_i}{\sum_i X_i/A_i},
\qquad
\bar Z = \bar A \sum_i \frac{X_i}{A_i}Z_i,
\]

\[
Y_e = \sum_i X_i\frac{Z_i}{A_i},
\qquad
m_c = \sum_i X_i\frac{W_i}{A_i}.
\]

For a generic charge moment \(M_g=\bar A\sum_i X_i g_i/A_i\),

\[
\frac{\partial M_g}{\partial X_j}
=
\frac{\bar A}{A_j S_X}(g_j-M_g).
\]

The active Skye number-fraction derivative is computed in the same helper.
For active set \(A_{\rm Skye}\),

\[
y_i=\frac{X_i/A_i}{\sum_{a\in A_{\rm Skye}}X_a/A_a}.
\]

The active-set derivative is used by both Skye ideal ions and Skye Coulomb
composition partials so the two paths share the same constrained basis.

# Row-Scoped EOS API

The internal EOS API has two modes:

| Mode | Purpose | Rows |
|---|---|---|
| full dxa | diagnostics, plotter, finite-difference comparison | all `num_eos_basic_results` |
| row-scoped dxa | production star solve and Brunt faces | requested rows only |

Production wrappers include:

| Wrapper | Use |
|---|---|
| `get_eos_star_dxa` | cell-centered hydro rows |
| `get_eos_star_dxa_with_moments` | same, but reuses star-side moments |
| `get_eos_brunt_dxa` | face/path Brunt rows |
| `get_eos_brunt_dxa_with_moments` | face/path Brunt rows with precomputed moments |

The row-scoped path still writes into legacy `nv x species` buffers for stable
indexing, but it avoids computing expensive component rows that were not
requested.  A future cleanup can carry smaller row-local buffers internally.

Legacy `eosDT_get` callers, including split burn and net one-zone burn, stay on
`include_composition_partials=.false.` and do not enter the Skye/FreeEOS
composition-partial path.  They still need the two historical `d_dxa` rows, so
the public wrappers keep a local full-row scratch buffer and copy those rows
back without per-call heap allocation.

# Skye Thermodynamic Packing

Skye remains a Helmholtz-free-energy path.  For a composition derivative
\(F_j\), with \(T,\rho\) carried by `auto_diff_real_2var_order3`,

\[
S_j=-\frac{\partial F_j}{\partial T},
\qquad
P_{{\rm gas},j}=\rho^2\frac{\partial F_j}{\partial \rho},
\qquad
E_j=F_j+T S_j.
\]

Then

\[
\frac{\partial\ln E}{\partial X_j}=\frac{E_j}{E},
\qquad
\frac{\partial\ln P_{\rm gas}}{\partial X_j}
=\frac{P_{{\rm gas},j}}{P_{\rm gas}}.
\]

`eos/private/skye_thermodynamics.f90:pack_composition_partials` computes only
the requested rows.  Derived rows such as `grad_ad`, `gamma1`, and heat
capacity rows are available when requested but are not paid for in the normal
hydro-row path.

# Skye Ideal Ion Path

For

\[
F_{\rm ion}=\frac{k_B T}{m_u\bar A}\Phi,
\]

\[
\Phi = \sum_i y_i
\left[
\ln\left(\frac{y_i n}{n_{Q,i}}\right)-1
\right],
\qquad
n=\frac{\rho}{m_u\bar A},
\]

the implemented derivative reuses the active number-fraction helper:

\[
\Phi_j =
\sum_i y_{i,j}\ln\left(\frac{y_i n}{n_{Q,i}}\right)
-\frac{\bar A_j}{\bar A},
\]

\[
F_{{\rm ion},j}=
\frac{k_B T}{m_u}
\left(
\frac{\Phi_j}{\bar A}
-\frac{\Phi\bar A_j}{\bar A^2}
\right).
\]

The current speed refactor stores the repeated
\(\ln(y_i n/n_{Q,i})\) terms once per call and avoids rebuilding the full
active-fraction Jacobian for every isotope column.

# Skye Electron Path

Skye uses the HELM electron helper for the electron/positron free energy.  With

\[
D=Y_e\rho,
\qquad
F_e=Y_e f(T,D),
\]

the composition derivative of each \(T,\rho\) coefficient is

\[
\frac{\partial}{\partial X_j}
\left[
\frac{\partial^{a+b}F_e}{\partial T^a\partial\rho^b}
\right]
=
Y_{e,j}
\left[
(b+1)Y_e^b f_{a,b}
+ \rho Y_e^{b+1} f_{a,b+1}
\right].
\]

The Skye wrapper maps this through \(Y_{e,j}\) and fills electron/free-electron
rows only when requested.

# Skye Coulomb Path

## Full Companion Path

The full companion path exists for general composition partials, phase rows,
latent-heat rows, and phase-transition cells.  It uses
`skye_composition_ad_real`, which carries

\[
(u, u_X),
\]

where each entry is still a `auto_diff_real_2var_order3` in \(T,\rho\).  This
lets the existing Coulomb formulas be copied closely while carrying a single
composition perturbation.

It covers:

- liquid and solid OCP free energies;
- liquid and solid screening terms;
- liquid and solid mixing corrections;
- electron exchange-correlation;
- mixing entropy;
- Gamma-limit extrapolation;
- hard branch and soft phase transition.

## Production Hydro-Row Basis

For ordinary hydro rows away from the phase transition, the current production
path uses the basis

\[
\{Y_e,\bar A,y_1,\ldots,y_N\}.
\]

Network isotope columns are assembled by

\[
F_{C,j}
=
Y_{e,j}F_{C,Y_e}
+\bar A_jF_{C,\bar A}
+\sum_i y_{i,j}F_{C,y_i}.
\]

This reduces the expensive Coulomb work from one companion call per network
isotope to a small basis.

## Off-Transition Hard-Branch Reuse

When the base Skye EOS returns exactly liquid or exactly solid and phase/latent
composition rows are not requested, the hard-branch value is

\[
F_C = F_{\rm branch}(T,\rho,Y_e,\bar A,\{y_i\}).
\]

In this branch:

\[
x_{\rm nefer}=N_A Y_e\rho,
\]

so the \(Y_e\) basis derivative is exactly

\[
F_{C,Y_e}
=
\frac{\partial F_C}{\partial \rho}\frac{\rho}{Y_e}.
\]

The \(\bar A\) basis derivative enters the Coulomb correction through the final
unit conversion \(kT/(\bar A m_u)\), so

\[
F_{C,\bar A}=-\frac{F_C}{\bar A}.
\]

These are algebraic reuses of the base Skye AD result; they are not physics
approximations.  The full companion path remains active for phase-transition
cells and phase/latent rows.

## Batched Active-Number-Fraction Path

The active-number-fraction basis uses `nonideal_corrections_dya`.  For a
single selected branch, each OCP leaf value

\[
f_i(T,\rho)
\]

is independent of the active number-fraction perturbation.  Therefore the
branch sum

\[
F_{\rm OCP}=\sum_i y_i f_i
\]

has

\[
\frac{\partial F_{\rm OCP}}{\partial y_j}=f_j.
\]

The current implementation caches the selected branch leaf values from the base
Skye Coulomb evaluation and reuses them in the hard-branch batched `dYA` path.
Near the phase transition it falls back to the full evaluation.

# FreeEOS / OPAL / SCVH Table Components

FreeEOS and OPAL/SCVH naturally provide a compact table basis in \(X,Z\).
The current path keeps this compact pair through the component call and expands
to isotope columns with the shared constrained helper.

This avoids doing the table work once per isotope.  Timing from the 500-model
benchmark showed FreeEOS dxa cost around 11.4 thread seconds, with table lookup
around 2.6 thread seconds and `X,Z` expansion negligible.

# eosDT Blends

The current production policy is fixed-weight composition rows through eosDT
selector blends:

\[
R_j = \alpha R_{1,j} + (1-\alpha)R_{2,j}.
\]

The selector composition term

\[
\alpha_j(R_1-R_2)
\]

is intentionally not used for production composition rows.  This avoids adding
composition derivatives of numerical EOS selector functions to the physical
Jacobian.  Some selector derivative plumbing still exists for diagnostics and
future designs.

Long-term thermodynamic endpoint: if component EOSs can expose consistent
potential-level information, the cleaner blend would blend a Helmholtz free
energy first and derive \(P,E,S,\chi\), and composition rows from the single
blended potential.

# Brunt / Ledoux Path

The implicit Brunt path lives in `star/private/implicit_brunt.f90`.

With `implicit_diffusion_flag` and `use_Ledoux_criterion`, the stored Brunt
value defaults to the same two-composition finite pressure difference as the
ordinary Brunt path.  This keeps the sign seen by thermohaline activation on
the same value definition with and without implicit diffusion.  Setting
`implicit_diffusion_use_brunt_finite_difference_value = .false.` restores the
earlier stored value from the linearized EOS-partial contraction.  The implicit
Jacobian uses a linearized pressure-composition coefficient from face EOS
partials in both modes.

Schematically, the stored pressure-composition contribution is

\[
\Delta\ln P_X^{\rm val}
=
\ln P_{\rm eos}(\rho_f,T_f,X_k)
-
\ln P_{\rm eos}(\rho_f,T_f,X_{k-1}).
\]

The AD derivative carrier is the linearized contraction
\[
\sum_i \chi_{X_i,f}\Delta X_i,
\qquad
\chi_{X_i,f}
=
\left.
\frac{\partial\ln P_{\rm eos}}{\partial X_i}
\right|_{\rho_f,T_f,X_f}.
\]
The face EOS pressure-composition coefficient is treated as a current-iterate
coefficient.  The optional Gauss path changes this derivative coefficient,
not the stored Brunt value.  Second derivatives of EOS composition partials are
intentionally not part of this implementation slice.

When `use_Ledoux_criterion = .false.`, the implicit Brunt face EOS path is
skipped and the composition term is kept zero.  This is both physically
consistent with the Ledoux-off choice and avoids paying face EOS dxa costs.

# Implicit Diffusion Policy

The current split is

\[
D_{\rm mix} = D_{\rm mix}^{\rm implicit} + D_{\rm mix}^{\rm explicit}.
\]

In code, `D_mix(:)` is still the ordinary real MESA coefficient consumed by
existing mixing and sigma routines. `Dmix_implicit(:)` carries the AD
derivatives for the promoted component, and `Dmix_explicit(:)` is the real
semi-implicit coefficient held fixed through Newton iterations.

Implicit component:

- MLT/TDC convective mixing;
- semiconvection;
- thermohaline.

For a full `set_mixing_info` pass, the promoted set is selected from the
post-cleanup `mixing_type`.  For a solver iteration, the promoted set is
selected from current `mlt_mixing_type` after `set_mlt_vars`; matching promoted
`mixing_type` display flags are updated in the same path.

Explicit or semi-implicit component:

- overshoot;
- rotation;
- minimum mixing;
- `other_D_mix`;
- other post-processing or external hooks unless separately promoted.

The isotope residual uses the usual fixed-coefficient implicit diffusion form.
During a full `set_mixing_info` pass, `set_Dmix_components(s,.true.)` refreshes
both pieces directly: promoted MLT/TDC/thermohaline cells use `mlt_D_ad` in
`Dmix_implicit`, and `Dmix_explicit` stores the non-implicit part of the current
ordinary MESA coefficient.  If the MESA total is smaller than the promoted local
value, the promoted AD coefficient is scaled down so both stored components
remain non-negative. The total coefficient is then rebuilt as
`D_mix = Dmix_implicit%val + Dmix_explicit`, matching the full-pass MESA total.

During Newton iterations, `set_Dmix_components(s,.false.)` keeps
`Dmix_explicit` from the last full `set_mixing_info` pass fixed and refreshes
only `Dmix_implicit`.  The promoted coefficient value and derivatives come
from current `mlt_D_ad`, and the solver-iteration promoted set follows current
`mlt_mixing_type` for convection, semiconvection, and thermohaline.  The same
path updates matching promoted `mixing_type` display flags, clears stale
promoted flags when the current promoted coefficient is inactive, and refreshes
`D_mix_non_rotation`.  Nonlocal full-pass edits such as overshoot, minimum
mixing, user zeroing, and boundary cleanup remain in `Dmix_explicit` unless
separately promoted.  The solver refreshes scalar `sig(:)` values every
implicit Newton iteration; nonlinear coefficient Jacobian terms are gated:

The sigma value always comes from total `D_mix`.  The optional sigma
derivatives are built from `Dmix_implicit`, so the Jacobian only differentiates
the promoted implicit component.

For thermohaline, `turb:set_thermohaline` now passes the scalar AD Ledoux
composition term into `thermohaline:get_D_thermohaline_ad`.  Kippenhahn and
Traxler are evaluated in AD arithmetic.  Brown et al. keeps the real Newton
solve for the fitted `(l,\lambda)` values and differentiates the converged fit
by solving
\[
  J_{F,(l,\lambda)}
  \begin{bmatrix}\partial l\\\partial \lambda\end{bmatrix}
  =
  -
  \begin{bmatrix}\partial F_1\\\partial F_2\end{bmatrix}_{(l,\lambda)\ {\rm fixed}},
\]
then applies those derivatives to
\[
  \mathrm{Nu}_\mu =
  1 + \frac{49\lambda^2}{\tau l^2(\lambda+\tau l^2)} .
\]

| Control | Default | Meaning |
|---|---:|---|
| `implicit_diffusion_use_brunt_finite_difference_value` | true | use ordinary finite two-composition pressure difference for the stored implicit Brunt value |
| `implicit_diffusion_include_dsig_structure` | false | include structure derivatives of sigma |
| `implicit_diffusion_include_dsig_dxa` | false | include high-gain cross-species sigma derivatives |

The `dxa` sigma block is only meaningful with Ledoux on, because its current
source is the Ledoux composition contribution to Brunt.

# Runtime Timing Snapshot

Recent user timing for a 500-model, no-implicit-diffusion benchmark with EOS
composition partials enabled:

| Quantity | Value |
|---|---:|
| total wall time | 30.715 s |
| EOS wall time | 7.979 s |
| `thread_time_eos_dxa/threads` | 7.721 s |
| `thread_eos_dxa_FreeEOS` | 11.355 thread s |
| `thread_eos_dxa_Skye` | 47.303 thread s |
| `thread_skye_dxa_ideal` | 1.123 thread s |
| `thread_skye_dxa_coul` | 16.700 thread s |
| `thread_skye_dxa_pack` | 0.336 thread s |

The main gain came from replacing redundant Skye Coulomb companion calls for
`dYe` and `dabar` with algebraic reuse of the base AD result.  A later speed
refactor caches selected OCP leaf values for the hard-branch `dYA` path.

# Verification Status

Completed earlier:

- install completed after the HELM `mu` derivative fix;
- EOS composition plotter produced finite-difference comparison outputs;
- static sign review of the implicit diffusion residual matched the existing
  `set_dxdt_mix` convention.

Still pending:

- compile/install after the latest Skye Coulomb speed refactors;
- focused finite-difference checks for the hard liquid, hard solid,
  extrapolated, and phase-transition Coulomb states;
- Brunt/Dmix runtime validation on the target stellar cases;
- HELM/FreeEOS derived-row audit;
- pure-PC internal composition derivative path, if it becomes a goal.

# Open Follow-Ups

1. Keep the row-scoped API but eventually shrink internal `nv x species`
   scratch to row-sized blocks.
2. Add focused Skye Coulomb finite-difference tests for the new hard-branch
   reuse and OCP-cache path.
3. Decide whether selector composition derivatives should remain diagnostic
   only, or whether the project should move toward a potential-level blend.
4. Validate the implicit Brunt/Dmix path in the high-temperature burning cases
   that motivated the work.
