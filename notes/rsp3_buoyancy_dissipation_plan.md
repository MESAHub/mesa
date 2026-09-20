# RSP3 buoyancy dissipation without moving the moments

Date: 2026-09-20. Proposal only. No MESA source or case changes.
No compilation or model execution.

## Sources and scope

- Felix Ahlborn's thesis, printed pp. 70--71, equations (4.2)--(4.6).
  Local source: `/Users/owner/Downloads/Ahlborn_Felix_thesis.pdf`.
  Printed p. 70 was checked visually as well as through extracted text.
- [Kupka, Ahlborn and Weiss (2022)](https://arxiv.org/abs/2207.12296),
  section 3.6, equations (19)--(22).
- [Braun et al. (2026)](https://arxiv.org/html/2604.06151v1#S2),
  section 2, including the use of the length in other closures.
- [Kovacs et al. (2026)](https://arxiv.org/html/2601.04931v1#S3.S1),
  section 3.1. Its alpha_tau must not be confused with the existing
  constant radial beta control.

The user wants to retain cell w and face Pi/Phi. The limited first
candidate described here changes only the turbulent energy dissipation
length. This is an explicit restriction of the published common-length
prescription. It is not an implementation of all its effects on moments,
radiation, viscosity, transport and mixing.

## Existing harmonic length and number of controls

`harmonic_dissipation_length_beta` is a real control, not a logical.
It is already `1d0` in the Cepheid case. Denote its value by beta_r:

```math
 \Lambda_h^{-1}=(\alpha H_P)^{-1}+(\beta_r r)^{-1}.
```

`get_TDC_mixing_length_cell` and `get_TDC_mixing_length_face` already
obtain this length through `get_mlt_mixing_length`. Retain these public
paths for the unmodified lengths and EOS scale heights.

One new nonnegative real coefficient for c4 is sufficient. Zero disables
the extension, so a separate enable flag is unnecessary. A first version
would require positive `harmonic_dissipation_length_beta` and RSP3.
The control name is not yet selected. The physical normalization is

```math
 c_4=\frac{0.3}{1.92 C_D},\qquad
 C_D={\tt RSP2\_alfad}\frac83\sqrt{\frac23}.
```

For unit RSP2_alfad this gives 0.0717624, approximately 0.072.
The approximately 0.2 value quoted with a different kinetic dissipation
normalization is not the same default. If C_D is varied, either adjust
c4 consistently or explicitly treat it as an independent coefficient.
No automatic division by zero is permitted when RSP2_alfad is zero.

## Eliminate the apparent division by w

With beta_r=1 the published relations are

```math
 \beta_s=\left(1+c_4\Lambda_d N/w\right)^{-1},\qquad
 \Lambda_d^{-1}=(\alpha H_P)^{-1}+(\beta_s r)^{-1}.
```

For the existing arbitrary positive beta_r, a proposed extension is to
use beta_r*beta_s as the radial factor. This recovers the published
normalization at beta_r=1. Algebra gives

```math
 \frac{c_4 N}{\beta_r r}\Lambda_d^2
 +\frac{w}{\Lambda_h}\Lambda_d-w=0,
```

```math
 \frac{w}{\Lambda_d}
 =\frac12\left[
 \frac{w}{\Lambda_h}+
 \sqrt{\left(\frac{w}{\Lambda_h}\right)^2+
             \frac{4c_4Nw}{\beta_r r}}
 \right].
```

There is no iteration and no division by w in the last expression.
It is still essential to evaluate the actual residual products:

```math
 \frac{\mathrm{DAMP}}{w}=
 \frac{C_D}{2}\left[
 \frac{w^2}{\Lambda_h}+
 \sqrt{\frac{w^4}{\Lambda_h^2}+\frac{4c_4Nw^3}{\beta_r r}}
 \right],\qquad
 \mathrm{DAMP}=w\frac{\mathrm{DAMP}}{w}.
```

For fixed positive N and w approaching zero from above, DAMP/w is
proportional to w^(3/2), and DAMP to w^(5/2). Both have finite first
derivatives, which vanish at zero. At N=0 or c4=0 the ordinary harmonic
dissipation is recovered exactly. An exact-zero branch must return
the limiting value and derivatives without differentiating sqrt(0).
The treatment of negative Newton trial w must follow the existing
physical-domain handling; negative w cannot enter this radical.
No artificial positive w floor is needed for this dissipation product.

An 80-digit scalar check over w=1e-20 through 1e8 reproduced the original
harmonic and beta relations to relative error below 4e-80. This verifies
the algebra only, not the MESA derivatives or evolution.

## Start of step buoyancy frequency

A proposed first discretization freezes the stable frequency during
Newton iteration:

```math
 N=\sqrt{\max(N_{\rm start}^2,0)}.
```

Use the MESA EOS form of the signed Brunt frequency, including composition,
rather than the thesis's ideal-gas reduction. On faces this has the form

```math
 N^2=\frac{g^2\rho}{P}\frac{\chi_T}{\chi_\rho}
       \left[B-(\nabla_{\rm actual}-\nabla_{\rm ad})\right].
```

Using N^2>0 as the gate is a proposed Ledoux-stable extension. It agrees
with the published thermal gate for uniform composition but differs in
composition-stabilized superadiabatic layers. This distinction must be
documented. Do not combine a thermal-only gate with a composition-aware
N and assume it is continuous at the thermal boundary.

Do not substitute -Y directly. Dynamic gradL and the temperature equation
form affect that mapping. Reconstruct the actual start state gradient
consistently with `hydro_gradient_support` and the face EOS. For cell
dissipation, first map signed face N^2 to the cell, then select its
positive part. Derive this mapping and its boundary cases explicitly.

The superadiabatic limiter provides a precedent for holding a timescale
during Newton iteration, not a reusable stored frequency:

- `conv_time_scale` in `star_utils.f90` uses abs(N^2).
- `turb_info.f90` can replace that timescale with Lambda/mlt_vc.
- `tau_conv_start` is populated conditionally for the limiter.

Consequently use a separate derived cache for this coefficient, prepared
on the current mesh before the solve. It must not depend on profile
output flags or the superadiabatic limiter being enabled. Rebuild it
after load, restart, remeshing and any preparation that changes the
start state. A retry at unchanged start state should reproduce it.
No conserved interpolation or independent photo state is needed for a
cache that is always reconstructed before use.

Freezing N is a time lag. Current w and the dissipation product retain
their implicit derivatives. The lag avoids differentiating the stability
switch within the Newton solve, but is not a proof of pulsation stability
or unchanged growth rates. Check timestep convergence. LNA must explicitly
handle the response of N; using a plain-real cached N silently omits it.
The implementation must choose and document either the physical
continuous closure perturbation or a frozen-N approximation.

## Why changing the common length is a larger task

As w approaches zero at fixed positive N, the common shortened length
scales as sqrt(w). Thus w/Lambda scales as sqrt(w), while the radiative
coefficient 1/Lambda^2 diverges as 1/w. Phi may remain positive when w
is zero in our present implementation. These limits require additional
moment and radiation treatment if the shortened length is used there.
Freezing only N does not remove them.

The dissipation-only candidate leaves the existing lengths in Pi/Phi
relaxation, radiative cooling, Lt, Eq/Uq and mixing. It preserves the
current variable placement. A common-length implementation can also
retain that placement, but cannot be called the same small change.

## Implementation and validation checklist

- Add the c4 control through defaults, declarations and control IO.
- Add and prepare the derived stable-frequency cache with retry and mesh
  lifecycle checks. Preserve the off path.
- Apply the regular product in `hydro_rsp2:compute_D_div_w`.
  `compute_D` already derives DAMP from that product.
- Verify `compute_C` and both energy formulations use the same exchange,
  so turbulent dissipation heats the gas consistently.
- Update startup estimates and `star_LNA_turbulence_closures`, including
  the chosen treatment of the background buoyancy frequency.
- Check exact zero, small positive w, stable and unstable layers,
  composition gradients, face reconstruction and all temperature rows.
- Check both velocity flags, both energy forms, restart and both mesh
  paths. No Pi/Phi relocation is required.
- With user authorization, compile and perform controlled paired runs
  and a timestep comparison. Do not claim the broad tail or closure
  realizability problem is fixed before those tests.

Assessment: a contained but nontrivial change for dissipation alone.
The algebra is small; cache preparation, energy exchange and LNA define
most of the implementation and verification work. No implementation
has been authorized by this exploratory request.

## Model 1001 budget, 2026-09-20

The current 350-zone run has Lt=0, yet stable-layer turbulence persists.
At r=3.9784 Rsun the current code gives Eq=1.315384e-3,
Source=-5.020281e-4 and DAMP=1.102386e-8 erg/g/s. Their sum matches the
accepted energy increase. Buoyancy therefore already brakes the motion,
while the implemented eddy-viscous term supplies more energy locally.
Ordinary cascade dissipation is weak, but a missing beta length is not
the complete demonstrated cause. Test this budget in a paired comparison
before claiming that a shortened dissipation length fixes the tail.
Evidence: `output/review/rsp3_modes_profile_20260920/stable_layer_budget.json`.
