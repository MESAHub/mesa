# EOS composition partials plan

Status: planning note only.  No MESA compile, MESA run, or EOS table
regeneration has been done for this plan.

## Goal

Fill analytic composition derivatives for Skye and for the combined eosDT
result for an arbitrary network size.

The immediate practical target is the composition derivatives that `star`
actually uses:

```text
d_dxa(i_lnPgas, j) = d ln(Pgas) / d xa_j | rho,T
d_dxa(i_lnE,    j) = d ln(E)    / d xa_j | rho,T
```

The cleaner internal target is to compute `d_dxa(1:nv,1:species)` in the
EOS layer, then keep the public `num_eos_d_dxa_results = 2` export unless
we explicitly decide to audit and enlarge every `star` array that assumes
only two composition derivatives.

## Current code state

Legacy public EOS calls still use full internal composition derivative
scratch storage, but the hot wrappers keep it as local automatic scratch
rather than allocating it on every call:

```text
eos/public/eos_lib.f90:eosDT_get
  real(dp) :: d_dxa_eos(num_eos_basic_results, species)
  call Get_eosDT_Results(..., d_dxa_eos, ...)
  d_dxa(1:num_eos_d_dxa_results,:) =
      d_dxa_eos(1:num_eos_d_dxa_results,:)
```

but `eos/public/eos_def.f90` still sets:

```fortran
integer, parameter :: num_eos_d_dxa_results = 2
```

so only `i_lnPgas` and `i_lnE` are reliably returned to `star`.

Skye currently evaluates the EOS and then zeros all composition
derivatives:

```text
eos/private/skye.f90:Get_Skye_EOS_Results
  call skye_eos(..., d_dxa, ierr)
  d_dxa = 0
```

The ideal EOS and PC EOS do the same in `eos/private/ideal.f90` and
`eos/private/eospc_eval.f90`.

For tabulated eosDT/FreeEOS tables, `eos/private/eosdt_eval.f90:get1_for_eosdt`
already maps the table derivatives with respect to `X` and `Z` onto the
network species:

```fortran
case (1)
   d_dxa(:,i) = d_dX
case (2)
   d_dxa(:,i) = 0
case default
   d_dxa(:,i) = d_dZ
```

This is consistent with the way `star/private/star_solver.f90` forms
conserved-composition solver derivatives by subtracting a sink species:

```text
d/d(x_j - x_sink) = d/dxa_j - d/dxa_sink
```

The eosDT blend logic currently ignores composition derivatives of blend
weights:

```text
eos/private/eosdt_support.f90:Do_Blend
  d_dxa = alfa*d_dxa_1 + beta*d_dxa_2
```

The missing term is:

```text
d alfa/dxa_j * (res_1 - res_2)
```

with the same quintic smoothing derivative used for the `lnT` and `lnRho`
partials.

## Derivative and gauge convention

EOS composition partials are fixed-`rho,T` analytic derivatives returned
as a full species-basis vector:

```text
d_dxa(:,j) = d res(:)/d xa_j | rho,T
```

This means analytic fixed-`rho,T` EOS derivatives, not the finite-difference
patch in `star/private/hydro_eqns.f90`.  The whole point of this work is
to remove reliance on that patch for Skye/eosDT-covered cells.

The composition vector and every physical composition perturbation must
satisfy the mass-fraction constraint:

```text
sum_j xa_j = 1
sum_j delta xa_j = 0
```

Because of that constraint, the full derivative vector has a gauge freedom:
adding the same constant to every species derivative does not change any
physical constrained derivative.  The meaningful objects are differences
or contractions with a zero-sum composition perturbation.

`star/private/star_solver.f90` applies this constraint by subtracting a
sink species derivative.  If species `s` is the sink, the derivative used
by a solver perturbation with `delta xa_j = delta` and
`delta xa_s = -delta` is:

```text
D_j^sink f = partial f/partial xa_j - partial f/partial xa_s
```

This matches the partial-checking path in
`star/private/star_solver.f90:dfridr_func`, where the tested abundance
perturbation is paired with an opposite sink perturbation.  Therefore the
EOS layer should not internally renormalize the final derivative into a
particular constrained `sum(xa)=1` derivative; the projection belongs at
the solver interface.

Existing MESA composition helpers already pick a normalized-composition
gauge.  For example `chem/public/chem_lib.f90:composition_info` returns:

```text
dabar_dx(j) = abar*(A_j - abar)/(A_j*sumx)
dzbar_dx(j) = abar*(Z_j - zbar)/(A_j*sumx)
dmc_dx(j)   = W_j/A_j - mass_correction
```

These are not meant to support arbitrary off-simplex composition states.
They are a full species-basis representative whose constrained differences
are physically meaningful.  A useful sanity check for quantities that
depend only on normalized composition is:

```text
sum_j xa_j * dQ/dxa_j = 0
```

not `sum_j dQ/dxa_j = 1`.

If a future caller wants a different projection, for example distributing
the compensating mass fraction over a set of species with weights `w_l`,
then:

```text
D_j^w f = d f/dxa_j - sum_l w_l d f/dxa_l
sum_l w_l = 1
```

The EOS should provide a complete analytic species-basis vector so either
projection can be formed.

This distinction is important for the Brunt form that contracts EOS
composition partials with a composition gradient:

```text
B = -(1/chiT) * sum_j (partial lnP/partial xa_j) * dxa_j/dlnP
```

The physical composition gradient must satisfy:

```text
sum_j dxa_j/dlnP = 0
```

Because of that zero-sum gradient, adding the same constant to every
`partial lnP/partial xa_j` leaves `B` unchanged.  What matters is that the
EOS partials are real analytic derivatives and that the composition
gradient is taken in the constrained composition space.

When an EOS component internally normalizes an active subset of species,
the derivative must include that component's actual normalization, because
that is part of the value being differentiated.

## Mathematical reference

This section is the compact mathematical contract for the implementation.
The code sections below expand where each piece should live.

### Composition moments and gauge

Let the network have `N = species` species, with mass fractions
`X_i = xa_i`, mass numbers `A_i`, charges `Z_i`, and atomic weights `W_i`.
MESA assumes:

```math
S_X \equiv \sum_{i=1}^{N} X_i = 1.
```

Define:

```math
Y_i = \frac{X_i}{A_i}, \qquad
Y = \sum_i Y_i, \qquad
\bar A = \frac{S_X}{Y},
```

and for any charge moment `g_i`:

```math
M_g = \bar A \sum_i \frac{X_i}{A_i} g_i.
```

The common choices are:

```math
\bar Z = M_Z, \qquad
\overline{Z^2} = M_{Z^2}, \qquad
\overline{Z^{5/3}} = M_{Z^{5/3}}.
```

The normalized-composition derivatives used by
`chem/public/chem_lib.f90:composition_info` are:

```math
\frac{\partial \bar A}{\partial X_j}
   = \frac{\bar A}{A_jS_X}\left(A_j-\bar A\right),
```

```math
\frac{\partial M_g}{\partial X_j}
   = \frac{\bar A}{A_jS_X}\left(g_j-M_g\right),
```

and:

```math
m_c = \sum_i X_i\frac{W_i}{A_i}, \qquad
\frac{\partial m_c}{\partial X_j}
   = \frac{W_j}{A_j}-m_c.
```

These satisfy the gauge check:

```math
\sum_j X_j \frac{\partial Q}{\partial X_j}=0
```

for normalized-composition quantities `Q`.  For complete ionization:

```math
Y_e = \frac{\bar Z}{\bar A}
    = \sum_i X_i\frac{Z_i}{A_i},
```

so:

```math
\frac{\partial Y_e}{\partial X_j}
   = \frac{Z_j}{A_j} - Y_e.
```

If a caller uses a sink species `s`, the constrained derivative is:

```math
D_j^{\rm sink} Q
   = \frac{\partial Q}{\partial X_j}
   - \frac{\partial Q}{\partial X_s}.
```

More generally, with compensating weights `w_l`:

```math
D_j^{w} Q
   = \frac{\partial Q}{\partial X_j}
   - \sum_l w_l\frac{\partial Q}{\partial X_l},
\qquad
\sum_l w_l = 1.
```

Any contraction with a zero-sum composition perturbation is gauge
invariant:

```math
\sum_j \delta X_j = 0
\quad\Longrightarrow\quad
\sum_j \left(\frac{\partial Q}{\partial X_j}+C\right)\delta X_j
 =
\sum_j \frac{\partial Q}{\partial X_j}\delta X_j.
```

### Brunt term from EOS partials

The Brunt branch `origin/Brunt_B_from_eos_partials` wants the chain-rule
form:

```math
B_{\rm comp}
   = -\frac{1}{\chi_T}
      \sum_j \chi_{X_j}\frac{dX_j}{d\ln P}.
```

With EOS output in MESA's current basis:

```math
g_{j,f} =
\left.
\frac{\partial \ln P_{\rm gas}}{\partial X_j}
\right|_{\rho_f,T_f,X_f},
```

where the EOS is evaluated at the face state
`(rho_face,T_face,xa_face)`.  This is the quantity that must come from
the EOS analytically.  It should not be approximated by differencing two
EOS calls with neighboring cell compositions.

The total-pressure conversion at the same face is:

```math
P_{\rm eos} = P_{\rm gas}+P_{\rm rad},
```

```math
\chi_{X_j}
   =
   \left.\frac{\partial \ln P_{\rm eos}}{\partial X_j}\right|_{\rho,T}
   =
   \frac{P_{\rm gas}}{P_{\rm eos}}\,g_{j,f},
```

and:

```math
\chi_T
   =
   \left.\frac{\partial \ln P_{\rm eos}}{\partial \ln T}\right|_{\rho,X}
   =
   \frac{P_{\rm gas}}{P_{\rm eos}}
      \left.\frac{\partial \ln P_{\rm gas}}{\partial \ln T}\right|_{\rho,X}
   + \frac{4P_{\rm rad}}{P_{\rm eos}}.
```

The remaining finite difference is the spatial composition gradient, not
an EOS composition partial:

```math
\frac{dX_j}{d\ln P}
   \approx
   \frac{X_{j,k-1}-X_{j,k}}
        {\ln P_{{\rm eos},k-1}-\ln P_{{\rm eos},k}},
```

which should satisfy:

```math
\sum_j \frac{dX_j}{d\ln P}=0
```

up to interpolation and roundoff.  If mass corrections are enabled, the
same form as the existing MHM implementation adds:

```math
B_{\rm mass}
   =
   -\frac{\chi_\rho}{\chi_T}
    \frac{d\ln m_c}{d\ln P}.
```

Current MESA's default Brunt path does not use this analytic contraction.
`star/private/brunt.f90:do_brunt_B_MHM_form` computes the composition
term by differencing two EOS pressures at the same face `rho,T`, but with
the two neighboring cell compositions:

```math
\ln P_1 =
\ln P_{\rm eos}(\rho_f,T_f,X_k),
\qquad
\ln P_2 =
\ln P_{\rm eos}(\rho_f,T_f,X_{k-1}),
```

```math
B_{\rm MHM}
 =
 \frac{\ln P_1-\ln P_2}
      {(\ln P_{{\rm eos},k-1}-\ln P_{{\rm eos},k})\chi_{T,f}}.
```

Linearizing this pressure difference gives:

```math
\ln P_1-\ln P_2
\approx
\sum_j
\chi_{X_j,f}
\left(X_{j,k}-X_{j,k-1}\right),
```

so:

```math
B_{\rm MHM}
\approx
-\frac{1}{\chi_{T,f}}
\sum_j
\chi_{X_j,f}
\frac{X_{j,k-1}-X_{j,k}}
     {\ln P_{{\rm eos},k-1}-\ln P_{{\rm eos},k}}.
```

This is the same first-order chain-rule expression, but the EOS
composition contraction is approximated by a finite pressure difference.
The new route makes the approximation explicit and narrower: the EOS
returns analytic face partials `chi_X_j,f`, and only the spatial
composition gradient is reconstructed from neighboring cells.

`star/private/hydro_vars.f90` then sets:

```text
s%gradL_composition_term(k) = s%smoothed_brunt_B(k)
```

when `use_Ledoux_criterion` and `calculate_Brunt_B` are true, and
`star/private/turb_support.f90` uses:

```text
gradL = grada + gradL_composition_term
```

so replacing the MHM pressure-difference Brunt calculation improves the
composition term in the existing Ledoux gradient path.

### Skye free-energy chain rule

Skye already evaluates a specific Helmholtz free energy:

```math
F = F_{\rm ion}+F_e+F_C.
```

For each species, define an analytic composition derivative object:

```math
F_j =
\left.\frac{\partial F}{\partial X_j}\right|_{\rho,T}.
```

The existing Skye thermodynamic identities are:

```math
S = -\frac{\partial F}{\partial T},
\qquad
P_{\rm gas} = \rho^2\frac{\partial F}{\partial \rho},
\qquad
E = F+TS.
```

Therefore the composition derivatives are:

```math
S_j = -\frac{\partial F_j}{\partial T},
\qquad
P_{{\rm gas},j}
   = \rho^2\frac{\partial F_j}{\partial \rho},
\qquad
E_j = F_j+T S_j.
```

The logarithmic EOS outputs are:

```math
\frac{\partial \ln S}{\partial X_j} = \frac{S_j}{S},
\qquad
\frac{\partial \ln E}{\partial X_j} = \frac{E_j}{E},
\qquad
\frac{\partial \ln P_{\rm gas}}{\partial X_j}
   = \frac{P_{{\rm gas},j}}{P_{\rm gas}}.
```

Radiation has no composition derivative at fixed `rho,T`, but it remains
in the total `S`, `E`, and `P` denominators where Skye packs total
quantities.

### Skye ideal-ion term

For the active Skye species set `A`, Skye uses active-set normalized
number fractions:

```math
S_A = \sum_{i\in A} X_i,
\qquad
q_i = \frac{X_i}{S_A},
\qquad
N_A = \sum_{i\in A}\frac{q_i}{A_i},
\qquad
y_i = \frac{q_i/A_i}{N_A}.
```

For `i,j in A`:

```math
\frac{\partial q_i}{\partial X_j}
   = \frac{\delta_{ij}}{S_A}-\frac{X_i}{S_A^2},
```

```math
\frac{\partial N_A}{\partial X_j}
   = \sum_{i\in A}\frac{1}{A_i}
     \frac{\partial q_i}{\partial X_j},
```

```math
\frac{\partial y_i}{\partial X_j}
   =
   \frac{(\partial q_i/\partial X_j)N_A/A_i
         -(q_i/A_i)(\partial N_A/\partial X_j)}
        {N_A^2}.
```

The ideal-ion free energy is:

```math
F_{\rm ion}
   =
   \frac{k_B T}{m_u\bar A}\Phi,
```

```math
\Phi
   =
   \sum_i y_i
   \left[
      \ln\left(\frac{y_i n}{n_{Q,i}}\right)-1
   \right],
\qquad
n = \frac{\rho}{m_u\bar A},
\qquad
n_{Q,i} = n_Q W_i^{3/2}.
```

At fixed `rho,T`:

```math
\Phi_j
   =
   \sum_i y_{i,j}
      \ln\left(\frac{y_i n}{n_{Q,i}}\right)
   - \frac{\bar A_j}{\bar A},
```

and:

```math
F_{{\rm ion},j}
   =
   \frac{k_B T}{m_u}
   \left(
      \frac{\Phi_j}{\bar A}
      - \frac{\Phi \bar A_j}{\bar A^2}
   \right).
```

If ion offsets are enabled:

```math
F_{\rm off}
   =
   \frac{\rm eV}{m_u}
   \sum_i \tilde y_i I(Z_i),
\qquad
F_{{\rm off},j}
   =
   \frac{\rm eV}{m_u}
   \sum_i \tilde y_{i,j} I(Z_i).
```

Here `tilde y_i` is the full-network normalized number fraction used by
`eos/private/ion_offset.f90:compute_ion_offset`, not the active Skye
Coulomb subset.

### Skye electron term

`eos/private/skye_ideal.f90:compute_ideal_ele` evaluates the HELM
electron free energy as a function of:

```math
D = Y_e\rho,
\qquad
F_e = Y_e f(T,D).
```

Let:

```math
f_{a,b}
   =
   \frac{\partial^{a+b}f}{\partial T^a\partial D^b}.
```

Then the composition derivative of the `T^a rho^b` coefficient is:

```math
\frac{\partial}{\partial X_j}
\left[
   \frac{\partial^{a+b}F_e}{\partial T^a\partial \rho^b}
\right]
 =
Y_{e,j}
\left[
   (b+1)Y_e^b f_{a,b}
   + \rho Y_e^{b+1} f_{a,b+1}
\right].
```

For HELM quantities that are functions of `T` and `D` but are not
multiplied by `Y_e`, for example `etaele`:

```math
Q_j = \rho Y_{e,j}\frac{\partial Q}{\partial D}.
```

`compute_ideal_ele` also returns a HELM `xnefer`, but current
`eos/private/skye.f90:skye_eos` overwrites it with
`compute_xne(den,ytot1,zbar)` before packing.  Therefore Skye's current
packed free-electron result uses complete ionization:

```math
x_{\rm nefer}=N_A\rho Y_e,
\qquad
\ln{\rm free}_e
   =
   \ln\left(\frac{x_{\rm nefer}}{N_A\rho}\right),
```

so at fixed `rho,T`, away from the `max(1d-99,xnefer)` clamp:

```math
\frac{\partial \ln{\rm free}_e}{\partial X_j}
   =
   \frac{Y_{e,j}}{Y_e}
   =
   \frac{\bar Z_j}{\bar Z}-\frac{\bar A_j}{\bar A}.
```

If the clamp is active, use the derivative of the clamped value, i.e. zero.

### Skye Coulomb term

The Coulomb term can be treated as a chain-rule problem over the scalar
composition-dependent inputs used by `skye_coulomb*.f90`:

```math
F_C = G(T,\rho,u_1,u_2,\ldots,u_m).
```

For each species:

```math
F_{C,j}
   =
   \sum_a
   \frac{\partial G}{\partial u_a}
   \frac{\partial u_a}{\partial X_j}.
```

The `u_a` include active-set number fractions, charge moments, `Y_e`,
`xnefer`, Coulomb parameters, mixing entropy, and phase-blend scalars.
This is why the plan uses a Skye-local companion type carrying:

```math
(u,\;u_1,\ldots,u_N)
```

where each entry is still an `auto_diff_real_2var_order3` in `T` and
`rho`.  Branches from `min`, `max`, and phase selection use the
derivative of the branch selected by the value path, away from exact
branch boundaries.

### Derived Skye result packing

For total pressure `P`, gas pressure `P_g`, and total energy `E`:

```math
C_v = \frac{\partial E}{\partial T},
\qquad
\chi_T = \frac{T}{P}\frac{\partial P}{\partial T},
\qquad
\chi_\rho = \frac{\rho}{P}\frac{\partial P}{\partial \rho}.
```

Composition derivatives are:

```math
C_{v,j} = \frac{\partial E_j}{\partial T},
```

```math
\chi_{T,j}
   =
   \frac{T}{P}\frac{\partial P_j}{\partial T}
   - \chi_T\frac{P_j}{P},
```

```math
\chi_{\rho,j}
   =
   \frac{\rho}{P}\frac{\partial P_j}{\partial \rho}
   - \chi_\rho\frac{P_j}{P}.
```

Using Skye's current derived-quantity identities:

```math
\gamma_3 - 1
   =
   \frac{P\chi_T}{\rho T C_v},
```

```math
\gamma_{3,j}
   =
   (\gamma_3-1)
   \left(
      \frac{P_j}{P}
      + \frac{\chi_{T,j}}{\chi_T}
      - \frac{C_{v,j}}{C_v}
   \right),
```

```math
\gamma_1
   =
   \chi_\rho + (\gamma_3-1)\chi_T,
```

```math
\gamma_{1,j}
   =
   \chi_{\rho,j}
   + \gamma_{3,j}\chi_T
   + (\gamma_3-1)\chi_{T,j},
```

```math
\nabla_{{\rm ad},j}
   =
   \frac{\gamma_{3,j}\gamma_1
         -(\gamma_3-1)\gamma_{1,j}}
        {\gamma_1^2},
```

and:

```math
C_{p,j}
   =
   C_p
   \left(
      \frac{C_{v,j}}{C_v}
      + \frac{\gamma_{1,j}}{\gamma_1}
      - \frac{\chi_{\rho,j}}{\chi_\rho}
   \right).
```

### eosDT table and blend equations

For eosDT/FreeEOS tables parameterized by hydrogen fraction `X_H` and
metallicity `Z`:

```math
R = R(T,\rho,X_H,Z).
```

The current isotope-class mapping is:

```math
\frac{\partial R}{\partial X_j}
 =
 \frac{\partial R}{\partial X_H}H_j
 + \frac{\partial R}{\partial Z}M_j,
```

where:

```math
H_j =
\begin{cases}
1, & Z_j=1,\\
0, & \hbox{otherwise},
\end{cases}
\qquad
M_j =
\begin{cases}
1, & Z_j>2,\\
0, & \hbox{otherwise}.
\end{cases}
```

Helium isotopes get zero table derivative in this gauge because the table
uses `Y = 1-X_H-Z`.

For a two-EOS blend:

```math
R = \alpha R_1 + (1-\alpha)R_2.
```

The complete composition derivative is:

```math
R_j =
\alpha R_{1,j}
 +(1-\alpha)R_{2,j}
 +\alpha_j(R_1-R_2).
```

If:

```math
\alpha = h(a_0)=10a_0^3-15a_0^4+6a_0^5,
```

then:

```math
\alpha_j = h'(a_0)a_{0,j},
\qquad
h'(a_0)=30(a_0-1)^2a_0^2.
```

### Optional coupled implicit-diffusion contribution

This should not be introduced as a second composition equation.  MESA
already solves one hydro composition equation for every isotope in
`star/private/hydro_chem_eqns.f90:do1_chem_eqns`.  In the current names:

```math
\mathrm{dxdt\_actual}_{j,k}
 =
 \frac{X^{n+1}_{j,k}-X^{n}_{j,k}}{\Delta t},
```

```math
\mathrm{dxdt\_expected}_{j,k}
 =
 \mathrm{dxdt\_mix}_{j,k}
 +
 \mathrm{dxdt\_nuc}_{j,k},
```

and the residual is:

```math
s\%\mathrm{equ}(i,k)
 =
 R^X_{j,k}
 =
 \frac{
   \mathrm{dxdt\_expected}_{j,k}
   -
   \mathrm{dxdt\_actual}_{j,k}
 }{\mathrm{eqn\_scale}_{j,k}}
 =0,
\qquad
i=s\%\mathrm{nvar\_hydro}+j.
```

If microscopic diffusion is moved into the same Newton solve, the
composition equation should keep this residual and add a diffusion rate to
`dxdt_expected`:

```math
\mathrm{dxdt\_expected}_{j,k}
 =
 \mathrm{dxdt\_mix}_{j,k}
 +
 \mathrm{dxdt\_nuc}_{j,k}
 +
 \mathrm{dxdt\_diff}_{j,k}.
```

For a mixed explicit/split plus coupled implementation, split the
transport coefficient into two accumulators:

```math
D_{\rm mix}
 =
D_{\rm mix}^{\rm explicit}
+
D_{\rm mix}^{\rm implicit}.
```

Here "explicit" means explicit with respect to the main hydro Newton
solve: it may remain in the current operator-split workflow, or it may be
added as a frozen source with no new Jacobian terms.  The implicit part is
the part whose value and partial derivatives are included in the Newton
residual:

```math
\mathrm{dxdt\_expected}_{j,k}
 =
 \mathrm{dxdt\_mix}_{j,k}
 +
 \mathrm{dxdt\_nuc}_{j,k}
 +
 \mathrm{dxdt\_diff,explicit}_{j,k}
 +
 \mathrm{dxdt\_diff,implicit}_{j,k}.
```

If the explicit part is kept operator-split, then
`\mathrm{dxdt\_diff,explicit}` is not added to this hydro residual; it is
applied by the existing split routine.  If it is added as a frozen source,
its Jacobian contribution is zero except for whatever frozen-source
bookkeeping MESA chooses to expose.

The conservative face-flux form is then the definition of the new rate,
not a new per-isotope equation:

```math
\mathrm{dxdt\_diff,implicit}_{j,k}
 =
 -\frac{
   \Phi^{\rm implicit}_{j,k+1/2}
   -
   \Phi^{\rm implicit}_{j,k-1/2}
 }{\Delta m_k},
```

The mass-fraction constraint can be enforced by either solving `N-1`
species with a sink species or by adding:

```math
R^\Sigma_k = \sum_j X_{j,k}-1=0,
\qquad
\frac{\partial R^\Sigma_k}{\partial X_{j,k}}=1.
```

The EOS composition partials enter the diffusion Jacobian wherever the
implicit fluxes depend on pressure, density, temperature, electron
fraction, mean molecular weight, Brunt terms, or thermodynamic response
functions.
For example:

```math
\frac{\partial R^X_{j,k}}{\partial X_{\ell,m}}
 =
 \frac{1}{\mathrm{eqn\_scale}_{j,k}}
 \left[
   \frac{\partial \mathrm{dxdt\_mix}_{j,k}}
        {\partial X_{\ell,m}}
 + \frac{\partial \mathrm{dxdt\_nuc}_{j,k}}
        {\partial X_{\ell,m}}
 + \frac{\partial \mathrm{dxdt\_diff,implicit}_{j,k}}
        {\partial X_{\ell,m}}
 - \frac{\partial \mathrm{dxdt\_actual}_{j,k}}
        {\partial X_{\ell,m}}
 \right].
```

The existing code already has the `dxdt_mix`, `dxdt_nuc`, and
`dxdt_actual` pieces.  The new coupled-diffusion work would be the
`dxdt_diff,implicit` value and its Jacobian entries.

For fluxes that need a molecular-weight derivative:

```math
\mu = \frac{\bar A}{1+\bar Z},
```

the composition derivative is:

```math
\mu_j =
\frac{\bar A_j(1+\bar Z)-\bar A\bar Z_j}{(1+\bar Z)^2},
\qquad
\frac{\partial\ln\mu}{\partial X_j} = \frac{\mu_j}{\mu}.
```

## Source-level verification against MESA

This verification is by direct source inspection only.  No MESA compile,
MESA run, test suite, or EOS table regeneration has been done.

The composition moment equations match
`chem/public/chem_lib.f90:composition_info`:

```fortran
sumx = sum(x(1:num_isos))
abar = sumx / sum(y(1:num_isos))
zbar = sum(y(1:num_isos)*z(1:num_isos)) * abar
z2bar = sum(y(1:num_isos)*z(1:num_isos)*z(1:num_isos)) * abar
...
z53bar = sum(y(i)*chem_isos%Z53(chem_id(i))) * abar
...
dabar_dx(i) = abar*(a(i)-abar)/a(i)/sumx
dzbar_dx(i) = abar*(z(i)-zbar)/a(i)/sumx
dmc_dx(i) = w(i)/a(i) - mass_correction
```

The general moment derivative in the math section is the same formula
with `g_i = Z_i`, `Z_i^2`, or `chem_isos%Z53`.  The `z53bar` equation
uses the final overwritten value based on `chem_isos%Z53`, not the
temporary `Z^3` assignment earlier in the routine.

The Skye ideal-ion equations match
`eos/private/skye_ideal.f90:compute_F_ideal_ion`:

```fortran
n = den / (amu * abar)
nQ = pow(kT, 1.5d0) / sifac
nj = ya(j) * n
nQj = nQ * pow(weights(j), 1.5d0)
F_ideal_ion = F_ideal_ion + ya(j) * (log(nj / nQj) - 1d0)
F_ideal_ion = F_ideal_ion * kt / (amu * abar)
```

The active-set normalization equations match
`eos/private/skye.f90:skye_eos`: species below
`mass_fraction_limit_for_Skye` are dropped, the remaining mass fractions
are normalized, and `ya(j) = select_xa(j)/A(j)` is normalized to sum to
one before the ideal-ion and Coulomb calls.

The free-energy thermodynamic equations match
`eos/private/skye_thermodynamics.f90:thermodynamics_from_free_energy`:

```fortran
s = -differentiate_1(F)
p = pow2(den) * differentiate_2(F)
e = F + temp * s
```

The derived-quantity equations match
`eos/private/skye_thermodynamics.f90:compute_derived_quantities`:

```fortran
chit = differentiate_1(p) * temp / p
chid = differentiate_2(p) * dens / p
cv = differentiate_1(e)
gam3 = 1d0 + (p / dens) * chit / (temp * cv)
gam1 = chid + (gam3 - 1d0) * chit
nabad = (gam3 - 1d0) / gam1
cp = cv * gam1 / chid
```

The radiation and packed Skye scalar equations match
`eos/private/skye_thermodynamics.f90:pack_for_export`:

```fortran
prad = crad * pow4(temp) / 3d0
erad = crad * pow4(temp) / dens
srad = 4d0 * crad * pow3(temp) / (3d0 * dens)
p = prad + pgas
e = erad + egas
s = srad + sgas
mu = abar / (1d0 + zbar)
lnfree_e = log(max(1d-99, xnefer)/(avo*dens))
```

The electron free-energy derivative formula matches
`eos/private/skye_ideal.f90:compute_ideal_ele`, which stores:

```fortran
F%val = free * ye
F%d1val1 = df_t * ye
F%d1val2 = df_d * pow2(ye)
F%d2val1 = df_tt * ye
F%d1val1_d1val2 = df_dt * pow2(ye)
F%d2val2 = df_dd * pow3(ye)
F%d3val2 = df_ddd * pow4(ye)
```

That is exactly:

```math
\frac{\partial^{a+b}F_e}{\partial T^a\partial\rho^b}
  = Y_e^{b+1} f_{a,b}.
```

The packed `lnfree_e` equation uses complete-ionization `xnefer`, because
`eos/private/skye.f90:skye_eos` calls `compute_ideal_ele` and then
overwrites `xnefer` with:

```fortran
xnefer = compute_xne(den, ytot1, zbar)
```

where `compute_xne` returns `avo*den*zbar/abar = avo*den*Y_e`.

The eosDT blend equation matches
`eos/private/eosdt_support.f90:Do_Blend` for `lnT` and `lnRho`
derivatives, and identifies the current missing composition term:

```fortran
res(j) = alfa*res_1(j) + beta*res_2(j)
dlnd(j) = alfa*d_dlnd_1(j) + beta*d_dlnd_2(j) +
          d_alfa_dlnd*res_1(j) + d_beta_dlnd*res_2(j)
d_dxa(j,k) = alfa*d_dxa_1(j,k) + beta*d_dxa_2(j,k)
```

Since `d_beta = -d_alfa`, the missing composition term is:

```math
\alpha_k(R_1-R_2).
```

The eosDT table mapping matches
`eos/private/eosdt_eval.f90:get1_for_eosdt`:

```fortran
case (1)
   d_dxa(:,i) = d_dX
case (2)
   d_dxa(:,i) = 0
case default
   d_dxa(:,i) = d_dZ
```

The gas-to-total pressure conversion for star-side partials matches
`star/private/micro.f90:store_eos_for_cell`:

```fortran
s% chiT_for_partials(k) =
   (s% Pgas(k)*d_dlnT(i_lnPgas) + 4d0*s% Prad(k))/s% Peos(k)
s% dlnPeos_dxa_for_partials(j,k) =
   s% Pgas(k)*d_dxa(i_lnPgas,j)/s% Peos(k)
```

The Brunt equations match the local branch
`origin/Brunt_B_from_eos_partials:star/private/brunt.f90` in structure:

```fortran
spatial_derivative_dX_dlnP =
   (s%xa(i,k-1) - s%xa(i,k)) / delta_lnP
B_term = B_term - d_eos_dxa(i_lnPgas,i)*spatial_derivative_dX_dlnP
s% brunt_B(k) = B_term / chiT
```

The note's total-pressure version is the corrected form needed before
moving that branch forward, because the branch currently uses
`d_eos_dxa(i_lnPgas,:)` and `d_eos_dlnT(i_lnPgas)` directly.

## Composition helper quantities

Add or reuse a helper that returns composition moments and their
derivatives for every network species.  The existing
`chem/public/chem_lib.f90:composition_info` already returns
`dabar_dx`, `dzbar_dx`, and `dmc_dx`.  Skye also needs derivatives of
`ye`, `z2bar`, `z53bar`, and active-set normalized number fractions.

For species `j`, define:

```text
A_j = chem_isos%Z_plus_N(chem_id(j))
Z_j = chem_isos%Z(chem_id(j))
Y_j = xa_j/A_j
Y   = sum_i Y_i
abar = 1/Y                         assuming sum(xa)=1
zbar = sum_i Y_i Z_i / Y
ye   = zbar/abar = sum_i xa_i Z_i/A_i
```

The existing MESA derivatives are:

```text
d abar/dxa_j = abar*(A_j - abar)/(A_j*sumx)
d zbar/dxa_j = abar*(Z_j - zbar)/(A_j*sumx)
```

and the corresponding complete-ionization electron fraction derivative is:

```text
d ye/dxa_j = d(zbar/abar)/dxa_j
            = dzbar_j/abar - zbar*dabar_j/abar^2
            = Z_j/A_j - ye                 for sumx = 1 in this gauge
```

For active-set normalization in Skye, let `A` be the set of species with
`xa > mass_fraction_limit`:

```text
S_A = sum_{i in A} xa_i
q_i = xa_i/S_A
N_A = sum_{i in A} q_i/A_i
y_i = (q_i/A_i)/N_A
```

For `i,j in A`:

```text
dq_i/dxa_j = delta_ij/S_A - xa_i/S_A^2
dN_A/dxa_j = sum_i (dq_i/dxa_j)/A_i
dy_i/dxa_j = ((dq_i/dxa_j)/A_i * N_A - (q_i/A_i)*dN_A/dxa_j)/N_A^2
```

For `j not in A`, `dq_i/dxa_j = dy_i/dxa_j = 0` as long as the active
set does not change.  The derivative is discontinuous at the
`mass_fraction_limit` boundary; tests should avoid points exactly on that
boundary.

## Skye free-energy derivative strategy

Keep Skye's current `auto_diff_real_2var_order3` path for `T` and `rho`.
Do not try to introduce an arbitrary-size global autodiff type.

Instead compute, for each species `j`, an `auto_diff_real_2var_order3`
object:

```text
F_x(j) = dF/dxa_j | rho,T
```

The value field is `dF/dxa_j`, and its `T`/`rho` derivative fields are
the mixed derivatives needed by the existing thermodynamic machinery.  In
other words, `differentiate_1(F_x(j))` is
`d/dT(dF/dxa_j)`, and `differentiate_2(F_x(j))` is
`d/drho(dF/dxa_j)`.

Once `F_x(j)` is known, the base thermodynamic composition partials are
straightforward:

```text
s_x(j) = - dF_x(j)/dT
p_x(j) = rho^2 * dF_x(j)/drho
e_x(j) = F_x(j) + T*s_x(j)

d lnS/dxa_j    = s_x(j)/S
d lnE/dxa_j    = e_x(j)/E
d lnPgas/dxa_j = p_x(j)/Pgas
```

These equations should live in a Skye helper near
`eos/private/skye_thermodynamics.f90:pack_for_export`, for example:

```text
pack_composition_partials(...)
```

That keeps the result-packing logic close to the existing value packing.

## Skye ideal-ion contribution

For the ideal-ion part:

```text
F_ion = (kT/(amu*abar)) * Phi
Phi = sum_i y_i * (ln(y_i*n/nQ_i) - 1)
n = rho/(amu*abar)
nQ_i = nQ * W_i^(3/2)
```

At fixed `rho,T`, for species `j`:

```text
dPhi_j = sum_i dy_i_j * ln(y_i*n/nQ_i) - dabar_j/abar
dF_ion_j = (kT/amu) * (dPhi_j/abar - Phi*dabar_j/abar^2)
```

The `auto_diff_real_2var_order3` fields of `F_ion_j` are obtained by
evaluating the same expression with the existing `temp` and `den`
autodiff values.  This preserves the mixed `T`/`rho` derivatives required
for `lnE`, `lnPgas`, and derived EOS quantities.

If `Skye_use_ion_offsets` is true, add the ionization offset derivative:

```text
offset = (ev2erg/amu) * sum_i yfull_i I(Z_i)
d offset_j = (ev2erg/amu) * sum_i dyfull_i_j I(Z_i)
```

where `I(Z)` is the ionization table in `eos/private/ion_offset.f90` and
`yfull_i` is the full-network normalized number fraction used by
`compute_ion_offset`, not the active Skye Coulomb subset.

## Skye electron contribution

`eos/private/skye_ideal.f90:compute_ideal_ele` already evaluates the
HELM electron free-energy table as a function of:

```text
din = ye*rho
F_ele = ye * f(T, din)
```

It already has most of the table derivatives needed for `T` and `rho`,
and it has commented-out `abar`/`zbar` derivative storage for `etaele`
and the HELM-table `xnefer`.  Current `skye_eos` overwrites `xnefer`
with the complete-ionization value from `compute_xne` before packing, so
`lnfree_e` currently needs only `dye_dxa`.

For a generic coefficient with `a` temperature derivatives and `b`
density derivatives:

```text
F_ele_(T^a rho^b) = ye^(b+1) * f_(T^a din^b)
```

Therefore the composition derivative is:

```text
d/dxa_j F_ele_(T^a rho^b)
 = dye_j * ((b+1)*ye^b*f_(T^a din^b)
            + rho*ye^(b+1)*f_(T^a din^(b+1)))
```

For the immediate `lnE`/`lnPgas` target, we need only the value,
`T` derivative, and `rho` derivative of `F_ele_j`.  For all internal
`nv` composition derivatives, fill enough second-order fields to compute
`Cv`, `chiT`, `chiRho`, `gamma1`, `gamma3`, `grad_ad`,
`dE_dRho`, `dS_dT`, and `dS_dRho`.

Implementation options:

1. Extend `compute_ideal_ele` to optionally return `F_ele_dxa(:)`,
   and `etaele_dxa(:)`.  Add `xnefer_dxa(:)` only if the current
   `compute_xne` overwrite is removed.
2. Or add a sibling helper that reuses the same local table derivatives
   and returns only the composition derivative objects.

The sibling helper is less invasive but duplicates some table-evaluation
code.  Extending `compute_ideal_ele` is cleaner if the patch remains
readable.

## Skye Coulomb contribution

The Coulomb correction depends on composition through:

```text
y_i, Z_i moments, abar, ye, xnefer, RS, GAME, Smix,
solid/liquid free energies, phase softmin, and kT/(abar*amu)
```

Do not hand-expand every fit in `skye_coulomb_liquid.f90` and
`skye_coulomb_solid.f90` into separate species loops.  That would be
large, fragile, and hard to review.

Use a small Skye-local companion type instead:

```fortran
type skye_xdiff
   type(auto_diff_real_2var_order3) :: v
   type(auto_diff_real_2var_order3), allocatable :: dx(:)
end type
```

This is not a new global arbitrary-variable autodiff system.  It is a
local chain-rule carrier for composition derivatives of Skye scalars.
Overload or explicitly implement only the operations Skye Coulomb uses:

```text
+, -, *, /, log, exp, sqrt, pow, min/max branches used away from kinks
```

Then mirror the Coulomb value path with `skye_xdiff` inputs for:

```text
RS, GAME, Zmean, Z2mean, Z52, Z53, Z321, Smix, kT
```

and real constants for species charges and weights.  The result is:

```text
F_coul_x(:), phase_x(:), latent_ddlnT_x(:), latent_ddlnRho_x(:)
```

The phase softmin derivative should follow the current branch behavior:
inside the smooth transition, differentiate the logistic expression; in
the far branches, the derivative of `phase` is zero because the value path
sets `phase = 0` or `phase = 1`.

## Packing all Skye `d_dxa` results

After forming:

```text
F_gas_x = F_ideal_ion_x + F_coul_x + F_ele_x
```

compute:

```text
s_x = -differentiate_1(F_gas_x)
p_x = rho^2*differentiate_2(F_gas_x)
e_x = F_gas_x + T*s_x
```

Radiation has no composition derivative at fixed `rho,T`.

For each species:

```text
d_dxa(i_lnS,j)    = s_x%val / s%val
d_dxa(i_lnE,j)    = e_x%val / e%val
d_dxa(i_lnPgas,j) = p_x%val / pgas%val
```

For derived quantities:

```text
cv_x     = differentiate_1(e_x)
chiT_x   = T*differentiate_1(p_x)/p - chiT*p_x/p
chiRho_x = rho*differentiate_2(p_x)/p - chiRho*p_x/p

gamma3_x = (gamma3 - 1) * (p_x/p + chiT_x/chiT - cv_x/cv)
gamma1_x = chiRho_x + gamma3_x*chiT + (gamma3 - 1)*chiT_x
grad_ad_x = (gamma3_x*gamma1 - (gamma3 - 1)*gamma1_x)/gamma1^2
cp_x = cp * (cv_x/cv + gamma1_x/gamma1 - chiRho_x/chiRho)
```

and:

```text
d_dxa(i_Cv,j)       = cv_x%val
d_dxa(i_Cp,j)       = cp_x%val
d_dxa(i_chiT,j)     = chiT_x%val
d_dxa(i_chiRho,j)   = chiRho_x%val
d_dxa(i_gamma1,j)   = gamma1_x%val
d_dxa(i_gamma3,j)   = gamma3_x%val
d_dxa(i_grad_ad,j)  = grad_ad_x%val
d_dxa(i_dE_dRho,j)  = differentiate_2(e_x)%val
d_dxa(i_dS_dT,j)    = differentiate_1(s_x)%val
d_dxa(i_dS_dRho,j)  = differentiate_2(s_x)%val
```

For composition-only scalar results:

```text
mu = abar/(1 + zbar)
dmu_j = (dabar_j*(1 + zbar) - abar*dzbar_j)/(1 + zbar)^2

lnfree_e = ln(zbar/abar)
dlnfree_e_j = dzbar_j/zbar - dabar_j/abar
```

Use the existing `max(1d-99, xnefer)` clamp semantics: if the value path
is clamped, set the derivative to zero for that clamped term.

## eosDT blend derivatives

Extend the eosDT blend plumbing to carry composition derivatives of the
blend fraction:

```fortran
real(dp), dimension(species) :: d_alfa_dxa
```

Update:

```text
eos/private/eosdt_eval.f90:combine_for_eosdt
eos/private/eosdt_support.f90:Do_Blend
```

The result derivative should become:

```text
dres_j/dxa_k =
   alfa*dres1_j/dxa_k
 + beta *dres2_j/dxa_k
 + dalfa_dxa_k*(res1_j - res2_j)
```

If the blend uses the existing quintic smoothing

```text
alfa = 10 a0^3 - 15 a0^4 + 6 a0^5
A = d alfa/d a0 = 30*(a0 - 1)^2*a0^2
```

then:

```text
dalfa_dxa_k = A * da0_dxa_k
```

matching the existing treatment of `d_alfa_dlnT` and `d_alfa_dlnd`.

## eosDT alpha derivatives by component

Add composition derivative versions of the alpha routines in stages:

1. Components with no composition-dependent alpha:
   set `d_alfa_dxa(:)=0`.

2. OPAL/SCVH alpha:
   `get_opal_scvh_alfa_and_partials` currently depends on `Z` and on
   `logT/logRho`.  Add `d_alfa_dZ` where the region logic blends in `Z`,
   then map:

   ```text
   d_alfa_dxa_j = d_alfa_dZ * dZ/dxa_j
   dZ/dxa_j = 1 for metals, 0 otherwise
   ```

   matching the current table mapping.

3. Skye alpha:
   `Get_Skye_alfa_simple` has no composition dependence because it uses
   fixed `rq` thresholds.  `Get_Skye_alfa` depends on `abar` and `zbar`
   through the polygon bounds.  Add a companion routine returning
   `d_alfa_dabar` and `d_alfa_dzbar`, then:

   ```text
   d_alfa_dxa_j =
      d_alfa_dabar*dabar_dxa_j + d_alfa_dzbar*dzbar_dxa_j
   ```

   At `max(...)`, `min(...)`, and polygon corner branch boundaries, use
   the derivative of the branch selected by the value path.

4. PC alpha:
   `Get_PC_alfa` depends on the Coulomb parameter through `zbar/abar`.
   Add:

   ```text
   d ln gamma_e/dxa_j = (1/3)*(dzbar_j/zbar - dabar_j/abar)
   ```

   and carry it through the same scalar blend function used for
   `logT/logRho`.

5. CMS and FreeEOS:
   audit whether the active alpha depends on `X`, `Z`, `abar`, or `zbar`.
   Add the same moment-based derivative only where the value path has real
   composition dependence.

## eosDT table composition derivatives

The tabulated eosDT and FreeEOS tables are not arbitrary-composition
tables.  Their analytic composition derivative is the derivative of their
current interpolation coordinates, not new physics for every isotope.

Keep `Get1_eosdt_Results` returning:

```text
dres/dX and dres/dZ
```

Then keep mapping to arbitrary network species by isotope class:

```text
H isotope:     d/dxa_j = d/dX
He isotope:    d/dxa_j = 0
metal isotope: d/dxa_j = d/dZ
```

This is consistent with the current `X`, `Y`, `Z` table parameterization
and with `star` sink-species subtraction.  A more physical isotope-level
table derivative is impossible from the existing `X,Z` tables alone.

## Star and implicit diffusion implications

`mu` is currently a dependent EOS/profile quantity, not a hydro solver
variable.  The relevant code path is:

```text
eos/private/skye_thermodynamics.f90:pack_for_export
  mu = abar/(1 + zbar)
  res(i_mu) = mu

star/private/micro.f90:set_micro_vars
  s% mu(k) = res(i_mu)

star_data/public/star_data_step_work.inc
  real(dp), pointer :: mu(:)
```

The solver variables are the structure variables plus species mass
fractions.  `star_data/public/star_data_step_input.inc` defines `xa` as
the composition variables and stores structure-variable indices such as
`i_lnd` and `i_lnT`; there is no `i_mu` solver-variable index there.

The current element diffusion solve is implicit inside
`star/private/diffusion.f90:do_solve_diffusion`, but it is operator-split
from the main hydro solve.  `star/private/element_diffusion.f90` calls it
with the current `abar`, `ye`, `free_e`, `T`, `rho`, gradients, and `xa`,
then the diffusion solver updates `xa`.

The coupled version should reuse the existing isotope residual in
`star/private/hydro_chem_eqns.f90`, not add a separate species equation.
In code names, the current residual is:

```text
dxdt_actual = s%xa_sub_xa_start(j,k)/s%dt
dxdt_expected = dxdt_mix + dxdt_nuc
s%equ(i,k) = (dxdt_expected - dxdt_actual)/eqn_scale
```

or, in math:

```math
R^X_{j,k}
 =
 \frac{
   \mathrm{dxdt\_expected}_{j,k}
   -
   \mathrm{dxdt\_actual}_{j,k}
 }{\mathrm{eqn\_scale}_{j,k}}
 =
0,
\qquad
\mathrm{dxdt\_actual}_{j,k}
 =
\frac{X^{n+1}_{j,k}-X^n_{j,k}}{\Delta t}.
```

The fully coupled microscopic-diffusion contribution should enter as:

```text
dxdt_expected = dxdt_mix + dxdt_nuc + dxdt_diff
```

with the conservative rate:

```math
\mathrm{dxdt\_diff}_{j,k}
 =
 -\frac{
   \Phi_{j,k+1/2}-\Phi_{j,k-1/2}
 }{\Delta m_k}.
```

The sign of `Phi` should follow the convention chosen in the diffusion
implementation, but the important structural point is that `dxdt_diff` is
a new expected rate in the existing composition equation.

A more flexible implementation should split the transport coefficient
according to how each part is solved:

```text
D_mix = D_mix_explicit + D_mix_implicit
```

or, in math:

```math
D_{\rm mix}
 =
D_{\rm mix}^{\rm explicit}
+
D_{\rm mix}^{\rm implicit}.
```

The default MESA behavior corresponds to `D_mix_implicit = 0` for the new
coupled microscopic-diffusion route and all existing transport continuing
through the current paths.  A source moved into `D_mix_implicit` supplies
both a rate and Jacobian contribution:

```math
\mathrm{dxdt\_diff,implicit}_{j,k}
 =
 -\frac{
   \Phi^{\rm implicit}_{j,k+1/2}
   -
   \Phi^{\rm implicit}_{j,k-1/2}
 }{\Delta m_k},
```

```math
\frac{\partial R^X_{j,k}}{\partial X_{\ell,m}}
\supset
\frac{1}{\mathrm{eqn\_scale}_{j,k}}
\frac{\partial
 \mathrm{dxdt\_diff,implicit}_{j,k}}
{\partial X_{\ell,m}}.
```

The `D_mix_explicit` part can still be handled by the existing split
routine or as a frozen source.  In either case it does not require the new
EOS-composition Jacobian terms.  This lets one transport source be coupled
without forcing all rotation, ordinary mixing, or user-supplied diffusion
hooks into the same Newton block.

For a future fully coupled implicit diffusion implementation, the minimal
solver formulation should still keep composition as the independent
variable set.  If a diffusion flux or stability term needs a `mu`
gradient, form it from the EOS/composition partials:

```text
mu = abar/(1 + zbar)
dmu/dxa_j = (dabar_j*(1 + zbar) - abar*dzbar_j)/(1 + zbar)^2
d ln mu/dxa_j = (1/mu) * dmu/dxa_j
D_j^sink ln mu = d ln mu/dxa_j - d ln mu/dxa_s
```

### Keeping rotation and other transport opt-in

The current transport workflow has several independently controlled
sources:

```text
star/private/hydro_chem_eqns.f90
  dxdt_expected = dxdt_mix + dxdt_nuc

star/private/mix_info.f90:update_rotation_mixing_info
  s% D_mix = s% D_mix_non_rotation + s% D_mix_rotation

star/private/element_diffusion.f90:do_element_diffusion
  operator-split microscopic diffusion when s% do_element_diffusion

star_data/public/star_data_def.inc
  other_D_mix and other diffusion hooks
```

The lowest-risk implementation path is to keep all existing controls and
add a separate opt-in route for coupled microscopic diffusion.  Conceptual
controls could look like:

```text
do_element_diffusion = .true.                 ! current operator-split path
use_coupled_element_diffusion = .false.       ! new default-off Newton path
use_coupled_diffusion_rotation_terms = .false.
use_coupled_diffusion_other_D_mix = .false.
```

Internally that means accumulating separate transport pieces:

```text
D_mix_explicit(:) = D_mix_non_rotation(:) + D_mix_rotation(:) + ...
D_mix_implicit(:) = coupled microscopic-diffusion or selected coupled terms
D_mix(:) = D_mix_explicit(:) + D_mix_implicit(:)
```

The exact accounting may need more than one array if the coupled term is
not representable as a simple scalar diffusion coefficient, but the split
is the important interface: explicit/split terms keep the old path;
implicit terms provide residual and Jacobian contributions.

The exact names should follow MESA control style when implemented.  The
behavior should be:

1. Default: no behavior change.  `do_element_diffusion` still calls the
   existing operator-split solver, and rotation/ordinary mixing still
   enter through `s%dxdt_mix`.

2. Coupled element diffusion on: skip the operator-split `xa` update for
   the part being coupled, move that part from `D_mix_explicit` to
   `D_mix_implicit`, compute `dxdt_diff,implicit` and its Jacobian, and
   add it to `dxdt_expected` in `hydro_chem_eqns`.

3. Rotation and ordinary mixing remain in `dxdt_mix` unless a separate
   control explicitly asks to couple additional transport terms.  This
   matters because rotation currently changes `D_mix` through
   `D_mix_rotation` and also has angular-momentum equations; pulling all
   of that into the new diffusion Jacobian at once would broaden the
   project substantially.

4. Any coupled transport source must conserve total mass fraction in the
   same basis used by the solver:

   ```math
   \sum_j \mathrm{dxdt\_diff}_{j,k}=0
   \quad\hbox{or}\quad
   \sum_j \Phi_{j,k+1/2}=0,
   ```

   up to the chosen sink/projection convention.

5. `other_D_mix`, `use_other_diffusion`, and
   `use_other_diffusion_coefficients` should keep their current
   operator-split behavior unless an explicit coupled hook is added that
   returns both transport rates and partial derivatives.

## Local branch context

Two local remote branches are relevant:

```text
origin/Brunt_B_from_eos_partials
origin/implicit_diffusion
```

`origin/Brunt_B_from_eos_partials` adds a control named
`use_eos_partials_for_Brunt` and a new Brunt path in
`star/private/brunt.f90`.  The intended EOS-partial form is:

```text
spatial_derivative_dX_dlnP = (s%xa(i,k-1) - s%xa(i,k))/delta_lnP
B_term = B_term - d_eos_dxa(i_lnPgas,i)*spatial_derivative_dX_dlnP
s%brunt_B(k) = B_term/chiT
```

The branch also contains a temporary finite-difference fallback whose tip
commit is named:

```text
add broken finite difference for eos partials
```

and the code comments say it can be removed when EOS fully provides
composition partials.  That confirms the right direction here: implement
analytic `d_eos_dxa(i_lnPgas,:)` in the EOS, then use those partials
directly in the Brunt chain-rule form.  Do not preserve or extend the
finite-difference fallback as production behavior.

The important distinction is that `d_eos_dxa(i_lnPgas,i)` should be the
EOS derivative at the face state returned by `get_eos`:

```text
get_eos(..., xa_face(:,k), rho_face, logRho_face, T_face, logT_face, ...)
```

The expression `(s%xa(i,k-1)-s%xa(i,k))/delta_lnP` is only the spatial
composition gradient that multiplies the EOS partial.  It is not a
replacement for the EOS partial.

One detail in the branch needs to be handled when moving it forward: the
EOS returns `lnPgas`, while the Brunt expression should use total EOS
pressure derivatives when radiation pressure is present.  At the face,
the star-side conversion should be:

```text
Pgas_face = exp(res(i_lnPgas))
Peos_face = Pgas_face + Prad_face
chiT_face = (Pgas_face*d_eos_dlnT(i_lnPgas) + 4*Prad_face)/Peos_face
chiX_j = (Pgas_face/Peos_face)*d_eos_dxa(i_lnPgas,j)

B = -(1/chiT_face) * sum_j chiX_j * dxa_j/dlnP
```

This is the same gas-to-total pressure conversion used by
`star/private/micro.f90` for `s%dlnPeos_dxa_for_partials`.

`origin/implicit_diffusion` expands `auto_diff_real_star_order1` from 27
to 33 derivative slots and adds fixed slots for `mu` and `H1`:

```text
i_mu_m1, i_mu_00, i_mu_p1
i_H1_m1, i_H1_00, i_H1_p1
```

That branch is not an arbitrary-network composition derivative mechanism.
A production arbitrary-network design should not rely on one hard-coded
composition slot.  It should use EOS analytic composition partial arrays,
then project them into whatever constrained composition basis the star
solver or diffusion discretization uses.

## Implementation phases

1. Add `notes/eos_composition_partials_plan.md` and keep it updated as
   the implementation evolves.

2. Add a small composition helper, probably in a new private EOS module:

   ```text
   eos/private/eos_composition_partials.f90
   ```

   It should return `dabar_dxa`, `dzbar_dxa`, `dye_dxa`, active-set maps,
   and active number-fraction derivatives.

3. Implement Skye `i_lnPgas` and `i_lnE` composition derivatives first.
   This is the minimal change that removes the most important need for
   `star/private/hydro_eqns.f90:fix_d_eos_dxa_partials`.

4. Add eosDT blend `d_alfa_dxa` support in `Do_Blend`, initially with
   zeros for all alpha derivatives.  This is a mechanical API change that
   should preserve existing behavior.

5. Add nonzero `d_alfa_dxa` for Skye and OPAL/SCVH blend weights.  These
   are the highest-priority blend derivatives for Skye/eosDT behavior.

6. Extend Skye packing to fill all internal `d_dxa(1:nv,:)` results.
   Keep public `num_eos_d_dxa_results = 2` unless there is a separate
   audited decision to expose more to `star`.

7. Revisit `origin/Brunt_B_from_eos_partials` after
   `d_eos_dxa(i_lnPgas,:)` is analytic for Skye/eosDT.  The branch's
   Brunt chain-rule form should then be able to use the EOS partials
   directly, with finite differences used only for verification.

8. If coupled implicit diffusion needs `mu` derivatives in `star`, avoid
   broadening `num_eos_d_dxa_results` until the call-site audit is done.
   Prefer either a narrow helper for `dmu_dxa` from composition moments or
   an explicitly expanded EOS derivative export with all affected arrays
   updated together.

9. Decide whether to add analytic derivatives for PC and ideal in the
   same pass.  The ideal EOS shares much of the Skye ideal-ion path, so it
   is a low-cost follow-up.  PC is more like Skye Coulomb and should be
   treated carefully.

10. After explicit permission, compile EOS and run targeted checks.

## Verification plan

No compile or run should be done without explicit permission.

When allowed, use progressively stronger checks:

1. Unit-level finite-difference checks through `eosDT_get_component` for
   `which_eos = i_eos_Skye`, away from active-set thresholds and blend
   boundaries.

2. Compare analytic `d lnE/dxa_j` and `d lnPgas/dxa_j` against centered
   finite differences for raw EOS partials:

   ```text
   xa_plus  = xa; xa_plus(j)  = xa(j) + eps
   xa_minus = xa; xa_minus(j) = xa(j) - eps
   dfdx ~= (f(xa_plus) - f(xa_minus))/(2 eps)
   ```

   Also check constrained finite differences, which are the derivatives
   used by `star` and preserve `sum_j xa_j = 1`:

   ```text
   xa_plus(j)  = xa(j) + eps; xa_plus(sink)  = xa(sink) - eps
   xa_minus(j) = xa(j) - eps; xa_minus(sink) = xa(sink) + eps
   dfdx_constrained ~= (f(xa_plus) - f(xa_minus))/(2 eps)
   d/d(x_j - x_sink) = d/dxa_j - d/dxa_sink
   ```

3. Repeat across representative regimes:

   ```text
   ideal-ion dominated
   electron-degenerate but weak Coulomb
   liquid Coulomb
   solid Coulomb
   smooth liquid-solid transition, not exactly at branch endpoints
   eosDT Skye blend
   eosDT OPAL/SCVH-to-HELM or OPAL/SCVH-to-Skye blend
   ```

4. Verify that old behavior is unchanged when Skye/ideal/PC fractions are
   zero.

5. After Skye analytic derivatives are trusted, test with
   `fix_d_eos_dxa_partials = .false.` in a controlled case, but only after
   explicit permission to run MESA models.

## Main risks

1. Active-set discontinuities from `mass_fraction_limit_for_Skye`.
   Derivatives are only meaningful while the active set is unchanged.

2. Blend-boundary kinks from `min`, `max`, polygon distance branches, and
   phase selection.  Use the derivative of the selected value branch and
   avoid exact boundaries in derivative tests.

3. Public API mismatch.  Internally we can fill `nv` derivatives, but
   changing `num_eos_d_dxa_results` beyond 2 requires a broad `star`
   allocation and call-site audit.

4. Coulomb derivative complexity.  The local companion type keeps this
   manageable; hand-expanding every Coulomb fit would be too invasive.

5. Table EOS limitations.  eosDT tables only know `X` and `Z`; isotope-
   resolved composition derivatives cannot be recovered from those tables.
