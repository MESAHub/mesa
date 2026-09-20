# RSP2 and RSP3 face turbulent energy plan

2026-09-20. Design and source audit, not an implemented layout change.
Audited on `EbF/star_lna`, after checkpoint `6d6bc2848`, with the current
uncommitted RSP2 and LNA corrections. No MESA compilation or stellar run
is part of this update.

## Decision and scope

The proposed RSP3 layout puts `w`, `Pi` and `Phi` on the same faces as
`Y_face` and luminosity. Keep `w=sqrt(e_t)` as the solver variable, with
`e_t` in erg/g. Moving its location does not require another equation or
another physical alpha. Gas density, temperature and composition stay in
cells. `v_flag` velocities stay on faces and `u_flag` velocities stay in
cells.

The preferred eventual architecture is one face energy implementation for
both RSP2 and RSP3, with different convection closures. Validate the RSP3
conversion first and compare a face RSP2 implementation against the current
cell RSP2 before adopting it. The current RSP2 calculation is a useful
reference. Its lack of independent Pi/Phi means that the demonstrated RSP3
averaging defect is not evidence that ordinary RSP2 needs the same repair.

This is a substantial discretization change. It includes energy storage,
pressure work, viscosity, turbulent transport, boundaries and saved state.
It is not a substitution of `w(k)` into routines that currently expect a
cell value. Keep TDC hydro separate. Reuse its established face scale-height
functions and verified identities, without merging its closure or changing
its time weights.

The existing [layout audit](rsp3_layout_audit.md) contains the competing
layouts, stencil checks and model-1001 evidence. This document is the
implementation plan. The [covariance correction](rsp3_covariance_closure.md)
and [buoyancy length proposal](rsp3_buoyancy_dissipation_plan.md) are separate
physical changes. Collocation does not implement either one or cure the
homogeneous closure failure.

## What the placement change addresses

Currently cell energy receives contributions from Pi on both bounding
faces. Face Pi/Phi then use energy reconstructed from adjacent cells. In
a uniform background an alternating Pi pattern cancels from cell buoyancy
production, although it still enters the face moment and gas heat equations.
The saved model has this cancellation in the oscillating stable layers.

With face energy, the local source at face k uses Pi(k) directly. Pi(k)
and Phi(k) use w(k)^2 directly. The gas still receives the difference of
the physical convective luminosities at its two bounding faces. This
removes this particular double averaging. It does not require Lc=0 when
Y<0: a time-dependent signed heat flux can be physical in stable layers.
Its spatial pattern, phase and amplitude must converge with resolution.

## 1. Storage, indexing and cell energy

Use MESA's inward indexing throughout. Gas cell k lies between outer face
k and inner face k+1. There are nz gas cells and nz+1 mathematical faces.
The existing solver has slots for faces 1:nz; the final boundary value
must have an explicit rule rather than an out-of-range solver entry.

| Quantity | Location in the proposed layout | Meaning |
| --- | --- | --- |
| w, w_start | Face | Square root of specific turbulent energy, cm/s |
| Pi | Face | Full velocity/entropy covariance, erg cm/(g K s) |
| Phi | Face | Full entropy variance, erg^2/(g^2 K^2) |
| Y_face, gradL, gradT | Face | Retain the existing temperature-row mapping |
| Lr, Lc, Lt, L | Gas face | Luminosities entering the gas energy balance |
| Turbulent transport between moment volumes | Gas cell center | Distinct from the gas-face Lt diagnostic |
| rho, T, e, composition | Cell | Existing gas state |
| Ptrb and turbulent energy used in gas work/storage | Cell | Derived from face energy consistently |

For the initial quadrature choose

```math
 e_{t,\mathrm{face},k}=w_k^2,\qquad
 e_{t,\mathrm{cell},k}=\frac{w_k^2+w_{k+1}^2}{2}.
```

The face mass for this quadrature is

```math
 \Delta m_{\mathrm{face},1}=\Delta m_1/2,\qquad
 \Delta m_{\mathrm{face},k}=(\Delta m_{k-1}+\Delta m_k)/2,
 \qquad
 \Delta m_{\mathrm{face},nz+1}=\Delta m_{nz}/2.
```

Including both boundary halves gives the exact storage identity

```math
 \sum_{k=1}^{nz}\Delta m_k e_{t,\mathrm{cell},k}
 =\sum_{k=1}^{nz+1}\Delta m_{\mathrm{face},k}w_k^2.
```

The specific cell energy is the average of squares, not the square of
averaged w. TDC uses the same quadrature with its different normalization:
`0.75*(mlt_vc(k)**2 + mlt_vc(k+1)**2)` equals this expression when
`mlt_vc=sqrt(2/3)*w`.

These are mass-lumped face degrees of freedom. The quadrature is a design
choice, not a claim that a point value is an exact volume average on every
unequal mesh. Its accuracy on sharply varying dm must be checked. Use its
same masses in conservation and remapping. Thermodynamic interpolation
weights need not be these quadrature weights.

**Inner boundary:** with prescribed w(nz+1)=0, the last cell energy is
w(nz)^2/2. With an explicitly selected zero-gradient envelope condition
w(nz+1)=w(nz), it is w(nz)^2 and the last stored face carries an extra
half-cell mass. These are different discretizations. Existing rotation
`remap1_face_average2` in `tdc_hydro_support` uses the latter end-volume
geometry. It cannot be copied unchanged for a fixed-zero inner moment.
`star_utils:set_dm_bar` also has RSP-specific boundary behavior; do not
infer all end weights from its name.

## 2. Residual equations

For an active face, write the energy residual in integrated specific-energy
units before applying numerical scaling:

```math
 R_{w,k}=w_k^2-w_{k,\mathrm{start}}^2
 +\Delta(P_{\rm trb}dV)_{\mathrm{face},k}
 +\Delta t\frac{L_{t,\mathrm{center},k-1}^{\theta}
                    -L_{t,\mathrm{center},k}^{\theta}}
                   {\Delta m_{\mathrm{face},k}}
 -\Delta t\left(\mathrm{SOURCE}_k-\mathrm{DAMP}_k
                         -\mathrm{DAMPR}_k+E_{q,k}\right)=0.
```

`Lt_center` here labels the proposed transport luminosity at a gas cell
center, not an existing array. The outer and inner boundary luminosities
replace the absent center terms at the domain ends. The pressure-work
increment is derived in section 5. It is not the old cell expression with
a face density substituted.

For RSP3, the proposed collocated energy source and existing decay are

```math
 \mathrm{SOURCE}_k=
 -\frac{1}{\rho_k}\left(\frac{\partial P}{\partial r}\right)_k
 \frac{\chi_{T,k}}{\chi_{\rho,k}c_{P,k}}\Pi_k,
 \qquad
 \mathrm{DAMP}_k=\texttt{RSP2\_alfad}\,
 \frac83\sqrt{\frac23}\frac{w_k^3}{\Lambda_k},
 \qquad \mathrm{DAMPR}_k=0.
```

All coefficients in this expression are face values. Preserve the pressure
gradient and EOS implementation in `rsp2_buoyancy_face`, including its AD
dependencies. Do not replace it by a hydrostatic gravity estimate. The
old `rsp2_moment_source` factors w_cell/w_face are no longer needed.

For a placement-only comparison, Pi and Phi retain their current physical
right-hand sides, evaluated with local w rather than reconstructed w:

```math
 R_{\Pi,k}=\Pi_k-\Pi_{k,\mathrm{start}}-\Delta t\left[
 \frac23 w_k^2\left(-\frac{\partial s}{\partial r}\right)_k
 -\frac{1}{\rho_k}\left(\frac{\partial P}{\partial r}\right)_k
       \frac{\chi_{T,k}}{\chi_{\rho,k}c_{P,k}}\Phi_k
 -\left(\texttt{RSP2\_alfa\_pi}\,6\sqrt{\frac23}\frac{w_k}{\Lambda_k}
       +\tau_{{\rm rad},k}^{-1}
       +\left(\frac{\partial u_r}{\partial r}\right)_k\right)\Pi_k
 \right]=0,
```

```math
 R_{\Phi,k}=\Phi_k-\Phi_{k,\mathrm{start}}-\Delta t\left[
 2\left(-\frac{\partial s}{\partial r}\right)_k\Pi_k
 -\left(\texttt{RSP2\_alfa\_phi}\,4\sqrt{\frac23}\frac{w_k}{\Lambda_k}
       +2\tau_{{\rm rad},k}^{-1}\right)\Phi_k\right]=0.
```

Here u_r means the mean radial velocity in either hydro formulation. Its
discrete strain follows the selected u/v implementation. These displayed
Pi/Phi equations intentionally describe the existing closure. They are
not asserted to preserve the covariance domain. The linked isotropic
closure proposal changes the Phi buoyancy coefficient in the Pi equation
to one third, ties Pi decay to the energy and Phi decay rates, and changes
the deformation term consistently. Test that physical correction separately
from placement, then test the combined formulation.

RSP3 keeps

```math
 L_{c,k}=4\pi r_k^2\rho_k T_k\Pi_k,\qquad
 R_{L,k}=L_{r,k}+L_{c,k}+L_{t,k}-L_k=0.
```

The current surface Y constraint and selected surface L/temperature row
remain explicit boundary rows. Moving w does not remove either row or
add a second total-luminosity variable.

## 3. Ordinary RSP2 on the face layout

RSP2 keeps its algebraic `PII` closure, not independent Pi/Phi. Evaluate
the existing `compute_PII_from_Y` at the same face as w. The candidate
face form of its current local source is

```math
 \mathrm{SOURCE}_k=(w_k+\texttt{RSP2\_source\_seed})
 \frac{P_k\chi_{T,k}}{\rho_k\chi_{\rho,k}c_{P,k}}
 \frac{\mathrm{PII}_k}{H_{P,k}}.
```

This relocates the factors in `compute_Source_div_w` and retains the
seed control's existing meaning. It is not automatically identical to
RSP3's pressure-gradient source, nor to the local equilibrium of its
two new equations. Those comparisons require matched closure coefficients.

The direct-face RSP2 luminosity is

```math
 L_{c,k}=4\pi r_k^2\frac{x_{\rm ALFAC}}{x_{\rm ALFAS}}
                (\rho T)_{\mathrm{face},k}\,w_k\,\mathrm{PII}_k.
```

Retain the existing reconstruction of rho*T in the placement comparison.
Replacing it by rho_face*T_face is a separate change on an unequal or
strongly stratified mesh. RSP2 radiative damping remains w_k^2/tau_rad,k;
RSP3 instead cools Pi/Phi. Preserve that difference and the alfar^2
normalization. Use `get_TDC_Hp_face` and `get_TDC_mixing_length_face` for
face closures, and the existing cell counterpart where a cell stress or
cell transport coefficient genuinely needs one.

Reasons to consider migrating RSP2 after validation are shared energy,
mesh and restart machinery and direct w coupling to PII/Lc. Reasons to
retain the current cell version until then are the changed spatial
quadrature, boundary behavior and pulsation growth. Fewer code paths alone
are not sufficient evidence to replace a working discretization.

## 4. Turbulent transport and the gas-face Lt

Moving energy to faces moves its nearest-neighbor diffusive flux to gas
cell centers. A candidate retaining the RSP2 alfat normalization is

```math
 L_{t,\mathrm{center},k}=
 -\texttt{RSP2\_alfat}\,(4\pi r_{\mathrm{center},k}^2)^2
   \rho_k^2\Lambda_{\mathrm{center},k}
   \frac{w_k+w_{k+1}}2
   \frac{w_k^2-w_{k+1}^2}{\Delta m_k}.
```

Outward luminosity is positive; indices increase inward. This uses an
arithmetic w for the mobility, with nonnegative w. Do not substitute the
cell energy gradient into the old face formula. Audit the exact center
radius and density quadrature against `compute_Lt`; changing those factors
is part of the convergence comparison.

The gas energy equation still needs Lt on its own bounding faces. For
interior faces choose

```math
 L_{t,k}=
 \frac{\Delta m_k L_{t,\mathrm{center},k-1}
       +\Delta m_{k-1}L_{t,\mathrm{center},k}}
      {\Delta m_{k-1}+\Delta m_k}.
```

At physical boundaries use the specified boundary transport luminosity.
With fixed masses and the section-1 half volumes, direct substitution gives

```math
 \left.\frac{d e_{t,\mathrm{cell},k}}{dt}\right|_{L_t}
 =\frac{L_{t,k+1}-L_{t,k}}{\Delta m_k}.
```

Thus turbulent transport cancels correctly when subtracting turbulent
storage from the combined gas/turbulent energy row. Compute the same
projected Lt in `compute_L_terms`, its caches, its start state and LNA.
Do not expose `Lt_center` as the existing profile `Lt` without changing
its documented meaning.

Keep the selected `L_theta_for_velocity_time_centering`. On a fixed mesh,
the linear projection commutes with this time weighting. After remesh,
old/start fluxes must be recomputed or consistently mapped before use.
This plan does not make total L fully implicit to suppress oscillations.

Nonzero alfat currently transports energy only. That can still violate
the joint moment bound by exporting energy while Pi/Phi remain. Location
alone cannot fix this. Future Pi/Phi transport belongs between the same
face control volumes, at gas cell centers. A common positive operator for
the three moments is a sufficient covariance-preserving spatial candidate;
arbitrary independent transport coefficients do not have that guarantee.
The moment transport fluxes have their respective moment units and are
not additional contributions to Lc. Keep this later physics extension
distinct from the placement test.

## 5. Turbulent pressure work and thermal energy

Use the projected cell energy in pressure:

```math
 P_{{\rm trb},k}=\texttt{RSP2\_alfap}\frac23\rho_k
                   \frac{w_k^2+w_{k+1}^2}{2}.
```

Preserve the existing pressure time weight, including rho_start and
w_start. For simple PdV work, the cell increment is

```math
 \Delta(P_{\rm trb}dV)_{\mathrm{cell},k}
 =P_{{\rm trb},k}^{\theta}
       (\rho_k^{-1}-\rho_{k,\mathrm{start}}^{-1}).
```

One globally conservative allocation splits this cell increment into the
contribution of each bounding face's energy. For j=k or k+1,

```math
 \Delta(P_{\rm trb}dV)_{k\leftarrow j}=
 \frac{\texttt{RSP2\_alfap}}3
 [\theta_P\rho_k w_j^2+(1-\theta_P)\rho_{k,\mathrm{start}}
                                      w_{j,\mathrm{start}}^2]
 (\rho_k^{-1}-\rho_{k,\mathrm{start}}^{-1}),
```

```math
 \Delta(P_{\rm trb}dV)_{\mathrm{face},j}=
 \frac{\sum_{k\text{ adjacent to }j}\Delta m_k
                   \Delta(P_{\rm trb}dV)_{k\leftarrow j}}
      {\Delta m_{\mathrm{face},j}}.
```

The sum of face work equals the sum of the original cell work, including
boundary contributions. No division by w is needed to form it.

This global identity is necessary but incomplete. Averaging face work
back into a cell does not generally reproduce that cell's original PdV
work. The difference can redistribute heat artificially unless it is
included in the consistent work-flux accounting. The implementation must
derive that local redistribution from the same face allocations, including
its compact Jacobian, or derive an alternative compatible pressure operator.
Do not declare conservation complete from a total-energy plot alone.

For both `eps_grav` and `dedt`, the combined energy row must use

```math
 \frac{d e_{t,\mathrm{cell},k}}{dt}
 =\frac{w_k^2+w_{k+1}^2-w_{k,\mathrm{start}}^2
                          -w_{k+1,\mathrm{start}}^2}{2\Delta t}.
```

Retain how `hydro_energy` distinguishes the forms:

| Energy/work selection | Required accounting |
| --- | --- |
| eps_grav | Gas thermodynamic pressure work stays in eps_grav; turbulent storage and turbulent pressure work remain explicit |
| dedt with simple PdV | Gas internal energy, turbulent storage and compatible pressure work |
| dedt with conservative mechanical work | Also retain the actual kinetic/potential storage and momentum-source work |

Check each with u_flag and v_flag. For RSP2 eps_grav, the current dispatch
uses the simple-work path even when the simple-work control is false;
preserve and document that dispatch. Gas pressure must not be counted
twice. The face conversion does not alter composition terms in eps_grav.
For changing cell masses, derive storage from changes of dm*e_cell rather
than applying the fixed-mass identities through a mass-transfer operation.

## 6. Eq and Uq for both velocity layouts

Use one stress in momentum and viscous heating. Preserve signs, geometry,
velocity time weights, mass corrections and the actual discrete mechanical
work. Keep the RSP2 routines separate from TDC's implementations.

Use TDC's eddy-viscosity discretization as the target precedent, while
keeping the implementations separate. Face w makes
that substantially simpler: both models then supply the velocity at the
same location, rather than requiring the current RSP2 cell/face ratios.

| Hydro choice | TDC spatial pattern to follow in RSP2/RSP3 |
| --- | --- |
| v_flag | Cell stress using the arithmetic average of its two face w values; face Uq from the stress difference; conservative allocation of Eq/w back to faces |
| u_flag | Face stress using local w; cell Uq from the two face stresses; local face Eq/w |

The bulk coefficient already has the same `(16/3)*pi*alpha_M` structure
in `tdc_hydro` and `hydro_rsp2`. Match `TDC_alpha_M` to `RSP2_alfam` when
comparing, and convert TDC `mlt_vc` to w by dividing by sqrt(2/3). This
does not require a new viscosity calibration parameter. Compare the stress,
force and heating separately on identical inputs before a stellar run.

Simplify the RSP2 viscosity routines in place: remove reconstructed w from
face stress, remove neighboring-w ratios used to divide by cell w, and use
the TDC-style allocation below. Keep the shared public scale-height calls.
Do not introduce a generic TDC/RSP2 viscosity wrapper during this change.

Matching the spatial pattern is not copying every TDC option. TDC can use
`mlt_vc_old` in momentum under
`TDC_alpha_M_use_explicit_mlt_vc_in_momentum_equation`; RSP2 uses its current
solver w. Preserve RSP2's current dependence and derivatives. Check masked
stress at the edge of a turbulent region, retained boundary force, inner
stress, pressure/radius time states and mass corrections explicitly. Those
details can differ even when interior formulas agree. Also distinguish the
native cell stress heating from the projected cell heating required by the
new turbulent storage below; their equality is an integrated identity,
not generally a pointwise one.

For v_flag, the stress remains in gas cells. A useful TDC precedent is
arithmetic face w in the cell viscosity, so that

```math
 E_{q,\mathrm{stress\ cell},k}
 =\frac{w_k+w_{k+1}}2\,(E_q/w)_{\mathrm{cell},k},
```

```math
 E_{q,\mathrm{face},j}=w_j
 \frac{\Delta m_{j-1}(E_q/w)_{\mathrm{cell},j-1}
       +\Delta m_j(E_q/w)_{\mathrm{cell},j}}
      {\Delta m_{j-1}+\Delta m_j}.
```

Use the corresponding one-sided expression at an end. Then

```math
 \sum_j\Delta m_{\mathrm{face},j}E_{q,\mathrm{face},j}
 =\sum_k\Delta m_k E_{q,\mathrm{stress\ cell},k}.
```

The cell `(Eq/w)` coefficient is formed analytically without dividing by
a possibly zero w. Its expression still uses the same cell stress and
new/time-centered velocity strains as Uq. This is different from replacing
the viscosity velocity by sqrt(e_cell): that alternative needs a new
allocation proof.

For u_flag, the stress and its heating already live on faces before the
current code averages heating into cells. Use local face w directly in
`compute_Chi_face` and retain the face force divergence in
`compute_Uq_dm_cell`. The old cell-division argument `k_div_w` cannot keep
its existing meaning. Check the inner-boundary stress separately, including
`TDC_include_inner_boundary_eddy_viscosity`.

In the combined gas/turbulent equation use
`Eq_cell = 0.5*(Eq_face(k)+Eq_face(k+1))`, matching the energy actually added
to the face equations. For v_flag this need not equal the native cell
stress dissipation point by point, although the integrated identity holds.
Retain the momentum work from the actual Uq in conservative dedt, including
the v_flag half-cell kinetic masses and optional mass corrections. Derive
and check the local stress-work redistribution as well as the global sum.

Do not assume the discrete Eq is nonnegative merely because continuum
viscosity dissipates energy: the present new-strain times midpoint-strain
product needs checking during a large timestep. This matters to the local
covariance proof. The zero-w limit and every AD shift must be checked.

## 7. Temperature rows and reconstruction

Preserve `hydro_gradient_support:get_rsp2_thermal_gradient` and the shared
coefficients of each temperature row. Y=gradT-gradL remains signed. Keep
the dynamical correction and composition correction in the same mapping,
with their current derivative and time-state treatment.

Coverage includes the gradient-difference form, the gradT form, dPrad/dm,
its opacity floor and radiation flux factor, dynamical gradL on/off,
composition gradients, and reconstructed face values on/off. No new
entropy-gradient smoothing is part of moving w. Reconstructed EOS and
scale heights still carry AD derivatives. Stored face w must not be
interpolated from adjacent w values a second time.

## 8. Zero moments, initialization and scaling

Moving w does not remove the zero derivative of w^2 at w=0. A face with
zero w, Pi and Phi, no incoming transport and no forcing may stay dormant.
A face with positive Phi and zero w is not generally dormant: Phi can
generate Pi and then turbulent energy. Do not force its three rows to zero.

The present divided energy residual, seed estimate and solver lower bound
were derived for cell w. Re-derive them for the face residual, including
transport and pressure work. Cancel known w factors analytically; do not
evaluate Pi/w or Eq/w at zero. An admissible local predictor must account
for the coupled Pi/Phi source. A fixed floor is not the derivation.

Initialization from MLT/TDC can use w_face=mlt_vc/sqrt(2/3) directly.
Initialization from an existing RSP2 model must remap energy, not average
w. If preserving initial Lc sets Pi=Lc/(4*pi*r^2*rho*T), check
Pi^2 <= (2/3)*w^2*Phi. That inequality constrains a seed; it is not a
unique equilibrium solution or a guarantee of equilibrium. Evaluate all
three initial residuals and report incompatible supplied moments.

Retain Pi/Phi's existing dimensional scales, recomputed from the correct
face state. Derive a face energy-rate scale from the neighboring existing
energy-rate scales with positive mass weights; freeze it for the solve.
The undivided row is scaled as Rw/(dt*energy_rate_scale_face). A divided
row requires the corresponding velocity factor, as in the current code.
Do not change L/Y normalization to compensate for the layout. Store and
check unscaled residuals so that convergence has the same physical meaning.

## 9. Boundary contract

Before enabling the new rows, write one table of allowed states for the
surface, center, truncated envelope, outer forced region and inner forced
region. Use that same mask for w, Pi/Phi, sources, fluxes, remesh and LNA.
Derive face limits from the cell-region controls; do not reuse an index
inequality merely because both arrays have length nz.

For each prescribed moment boundary specify its value, its quadrature
mass, the transport flux, the retained pressure/viscous boundary work and
where any imposed change of stored energy goes. A fixed-zero boundary
cannot silently discard a half-volume's remapped energy or receive
unaccounted viscous heating. Center symmetry and a truncated envelope are
different conditions. Preserve the existing supported boundary physics;
resolve the work allocation before imposing a blanket zero rule.

## 10. Remeshing in all three paths

Treat face moments as averages on their declared mass control volumes
for remapping. For positive overlap lengths use the same weights for all
three quantities:

```math
 (w^2,\Pi,\Phi)_{j,\mathrm{new}}=
 \sum_i \frac{\mathrm{overlap}_{ji}}
                    {\Delta m_{j,\mathrm{new}}}
                  (w^2,\Pi,\Phi)_{i,\mathrm{old}}.
```

Here overlap_ji is the mass shared by the new and old face intervals. Recover
w by the positive square root afterward. If the old covariance matrices
are positive semidefinite, this common convex remap preserves that property
and each moment's mass integral. It does not repair an inadmissible donor.
Although Pi/Phi are not conserved by physical source terms, their integrals
should not change simply because coordinates were remapped.

Start with piecewise-constant conservative overlap. A later higher-order
reconstruction must preserve positivity and the joint covariance, not just
monotonicity of three separate interpolants. Independent splines, or the
current interpolation of Pi/w followed by rescaling, do not give this
common-matrix guarantee. Y remains a signed point/interpolated face quantity;
it is not a conserved moment and must not use the energy overlap rule.

| Mesh path | Current code to replace or adapt |
| --- | --- |
| Ordinary adjustment | `mesh_adjust:do_etrb`, `do_RSP2_face_var`, and the thermal-energy update |
| Split/merge AMR | `adjust_mesh_split_merge` cell energy split/merge and Pi/w interpolation |
| Envelope remesh | `tdc_hydro_support:remesh_for_TDC_pulsations`, cell-energy overlap and moment reconstruction |

For split/merge, neighboring face control volumes change too. Save the
affected old dual intervals, perform a common overlap on the complete
affected patch, and include its boundary halves. Copying only the new
face at a split or deleting only one face at a merge is insufficient.
Consolidate the shared overlap operation only where it removes genuine
duplication across these paths; keep mesh geometry at its existing owners.

Check total thermal plus projected turbulent plus mechanical energy.
`mesh_adjust_get_T_from_E=.true.` must subtract the new projected turbulent
energy when recovering gas energy. The false path must use the same
definition and account for any intended thermal interpolation error; it
must not subtract face w(k)^2 as if it were a cell energy. The special
QHSE envelope remesher also changes hydrostatic structure: report its
thermal/gravitational adjustment separately from turbulent remap error.
Check artificial remap energy generation before blaming the next solve.

Rebuild Pi/Phi scales, Y/gradL-dependent quantities, pressure, Lt and mixing
after the new mesh and EOS are ready. No old-nz loop or stale start array
may be used. Changing mass domains through accretion or loss additionally
needs specified moment content of material entering/leaving the domain.

## 11. Saved models, photos and flag changes

The same `i_w` array cannot silently change from cell to face meaning.
The current photo layout records `RSP2_3equation_flag` from version 21;
the model file uses a separate RSP3 bit. Neither by itself identifies a
new face-w layout. Add versioned stored layout information to both formats
and restore it before interpreting xh. A stored flag is state metadata,
not a new user physics control. Inspect the current format version before
choosing a new value.

Keep exact restart of a new face photo distinct from conversion of an old
cell model/photo. Same-layout restart preserves moments and history.
Conversion maps energy onto the new control volumes and then rebuilds
dependent state. It cannot in general preserve every old cell energy,
every old Lc and every moment simultaneously. Match integrated energy,
report the local change and reject incompatible inputs rather than silently
clipping Pi or Phi. A conversion round trip can smooth a profile even
when its total energy is conserved.

When the selected layout changes, invalidate temporal extrapolation through
the existing generation machinery; initialize xh_start and w_start from
the converted state. The present RSP2/RSP3 flag-change routine preserves
gradT and redefines Y against the new gradL. Retain that ordering. If only
RSP3 is migrated initially, turning its flag off also requires a face-to-cell
energy conversion. Once both use faces, that extra conversion disappears.

No support for the obsolete Hp-equation layout is added. If old cell-w
restart conversion is deferred, fail with a clear layout diagnostic rather
than reading an old photo as face data. The known test photos are cell-w
snapshots and need this decision before any restart comparison.

## 12. LNA, Jacobian and diagnostics

LNA must linearize the selected spatial equations rather than reuse old
cell source/transport expressions with new variable names. Face energy
inertia contains 2*w_k*delta_w_k. The gas cell inertia contains

```math
 \delta e_{t,\mathrm{cell},k}
 =w_k\delta w_k+w_{k+1}\delta w_{k+1}.
```

Include the corresponding pressure/density inertia, work redistribution,
boundary rows and new Lt projection. Pi/Phi radiation and turnover losses
use the same face quantities as hydro. Eddy heating has no first-order
perturbation about zero mean strain, while viscous momentum damping remains;
do not apply that static simplification to a nonzero-strain background.

Audit `star_LNA_turbulence_closures` as well as `star_LNA_support`. Update
the term audit, work integrals, kinetic-energy quadrature and mode output
locations. Check the full-pencil residual with the actual new equations.
Do not use a residual tolerance of one to validate a conversion.

Do a structural stencil audit before coding. A face w row with a cell
stress coefficient can still be three-zone, but projected Lt in the gas
flux row can reach w(k-1:k+1), and velocity work/reconstruction can reach
farther. Existing `shift_p1`/`shift_m1` cannot retain an out-of-range AD slot.
Use MESA's established explicit partial accumulation where the physical
row requires it, or derive an equivalent compact row. Do not drop a
derivative or widen the physical stencil to conceal a missing partial.

Mixing uses v_conv=sqrt(2/3)*w_face and D=Lambda_face*v_conv/3 directly.
Keep the Y-based label criteria separate from actual mixing coefficients
and chemical mixing controls. Profile w and Pi/Phi have the same face
radius; expose cell turbulent energy separately where needed. Make the
meaning of the existing `etrb`, `Eq`, `SOURCE`, `DAMP`, `Lt` and history
integrals explicit. Remove all unintended second interpolation in outputs.

## 13. Source map

References are module/routine names so the plan survives line movement.
These are current code locations; proposed formulas above are not new APIs.

| File | Required audit/change |
| --- | --- |
| `star/private/hydro_rsp2.f90` | w row, moment RHS/source, PII/Lc, D/Dr, Lt, Eq/Uq, zero-state predictor, start state, moment remap |
| `star/private/auto_diff_support.f90` | w/etrb wrappers, explicit cell projection, start projection, convective velocity |
| `star/private/hydro_energy.f90` | `setup_d_turbulent_energy_dt`, all work forms, Eq and actual Uq work |
| `star/private/star_utils.f90` | `calc_Ptrb_ad_tw`, `cell_specific_total_energy`, integrated energy and boundary masses |
| `star/private/hydro_momentum.f90` | Shared stress forces, mass corrections and work consistency |
| `star/private/tdc_hydro.f90` | Existing public scale-height functions and face-energy viscosity precedent; preserve TDC behavior |
| `star/private/hydro_gradient_support.f90` | Preserve gradient mapping, reconstruction and derivatives |
| `star/private/hydro_eqns.f90` | Row dispatch, explicit partials and boundary temperature/L selection |
| `star/private/set_flags.f90` | MLT/RSP initialization and RSP2/RSP3 layout changes |
| `star/private/hydro_vars.f90` | Unpack xh, rebuild derived face/cell state |
| `star/private/solver_support.f90` | Scales, trial domain, zero-w handling and variable access |
| `star/private/alloc.f90`, `star_data/public/star_data_step_input.inc`, `star_data_step_work.inc` | Stored layout, arrays, current/old copies and any required cache |
| `star/private/photo_in.f90`, `photo_out.f90`, `read_model.f90`, `write_model.f90` | Versioned location and explicit old-layout conversion |
| `star/private/mesh_adjust.f90` | Common face-moment remap and both temperature recovery paths |
| `star/private/adjust_mesh_split_merge.f90` | Conservative affected-patch remap and total-energy subtraction |
| `star/private/tdc_hydro_support.f90` | Envelope remap, inner half-volume, QHSE and caches |
| `star/private/star_LNA_turbulence_closures.f90`, `star_LNA_support.f90` | Face residual/inertia, full partials, boundary rows, work and diagnostics |
| `star/private/profile_getval.f90`, `mix_info.f90`, `history.f90` | Location and normalization of outputs and energy integrals |

Keep the existing w wrappers as access to the stored unknown. Add a small
explicit cell-energy projection where shared consumers need it, with the
start-state counterpart. Do not globally redefine every `get_etrb` caller
without classifying whether it needs face or cell energy. Avoid new wrapper
layers that only rename existing calls.

## 14. Implementation sequence and acceptance

- [x] Audit current locations and the stable-layer averaging mechanism.
- [x] Specify candidate storage, transport and viscous projection identities.
- [x] Record the pressure-work redistribution and boundary requirements.
- [x] Check fixed-mass storage, Lt, Eq, pressure allocation and common overlap
      on unequal masses without MESA.
- [ ] Finish local work identities for u/v, both energy forms and all masks.
- [ ] Select the boundary contract and old-photo conversion behavior.
- [ ] Derive the face zero-state predictor and audit the complete AD stencil.
- [ ] Match the TDC viscosity spatial discretization in separate RSP2 routines,
      with coefficient, time-state, boundary and Eq/Uq identity checks.
- [ ] Implement RSP3 face storage, its closures, fluxes and gas projection.
- [ ] Wire every mesh, state, output and LNA path before stellar testing.
- [ ] Validate placement alone against the current closure on admissible states.
- [ ] Validate the separately derived local covariance correction.
- [ ] Compare RSP2 cell and face layouts with the same physical coefficients.
- [ ] Adopt a common face implementation only after both comparisons pass.
- [ ] Update the RSP3 and star_LNA manuscripts to the implemented equations.

Use a reviewable implementation checkpoint before replacing the current
layout. The existing checkpoint and current uncommitted corrections are
distinct states; do not discard either. No commit is made by this plan.
Prefer a branch comparison over adding a permanent experimental layout
control. If two layouts temporarily coexist, keep one clear stored layout
decision and remove transitional branches after validation.

Mathematical checks, without MESA:

1. Unequal masses: cell/face energy equality, boundary halves, Lt divergence
   equality, pressure allocation sum and Eq allocation sum.
2. Remap: identity mesh, constant state, narrow support, sign-changing Pi,
   split/merge, extreme mass ratios, total energy and covariance.
3. Smooth waves and alternating perturbations: retain direct face response;
   test accuracy rather than assuming that damping the alternating pattern
   proves convergence.
4. Zero-state and homogeneous source equations: domain preservation and
   finite derivatives, with/without radiation and energy transport.

MESA checks require the user's run/build authorization:

| Coverage | Evidence required |
| --- | --- |
| Both RSP2 and RSP3 | Startup, equilibrium residuals, periods, growth and work |
| Both u_flag and v_flag | Eq/Uq consistency and discrete mechanical energy |
| eps_grav and both dedt work forms | Local exchanges and integrated energy error |
| alfat=0 and nonzero | Turbulent flux cancellation; identify separate moment-domain failures |
| alfar=0 and 1 | Correct cooling normalization, stiff outer relaxation, no coupling hidden by cooling |
| alfap=0 and nonzero; alfam=0 and nonzero | Isolate storage, pressure and viscous terms |
| Every temperature row/reconstruction choice | AD partials including dynamic gradL and composition |
| Ordinary remesh, split/merge, envelope remesh | Before/after moment and total-energy budgets, both get_T_from_E settings |
| New photos and converted old inputs | Exact same-layout restart, explicit conversion error, retry rollback |
| Fixed 149/350/finer meshes | Tail sign changes, Pi/Phi phase, Lc, mode shape and mesh convergence |
| Several pulsation cycles | Accepted physical time, retries, timestep recovery, iterations and growth |
| Evolution beyond turnover times | Thermal diffusion limit, moving convective boundaries and composition gradients |

Do not select a replacement from a short run that merely avoids crashing.
Compare identical starting physical states, declared conversion errors and
elapsed stellar time. Keep luminosity time weighting fixed in the placement
comparison. Report cell/face conversion, local closure and added damping
as separate changes so their effects remain identifiable.

## Reproducible algebra checks

`notes/check_rsp2_face_w_plan.py` writes
`output/review/rsp3_face_w_plan_20260920/checks.json`. All assertions pass
for a mesh with maximum/minimum cell mass ratio about 7.4e5. Storage,
local Lt projection, Eq allocation, pressure allocation and remap integrals
agree within 7e-14 under the script's stated absolute/relative scaling.
Common remap keeps the covariance determinant positive. The script also
demonstrates a nonzero local pressure-work redistribution and smoothing
after a cell-to-face-to-cell conversion, despite conserved total energy.

The check includes both mathematical boundary faces. It does not impose
MESA's forced-zero masks, test AD partials, evolve the closure, validate
the full mechanical-work equations or run a star. Those checks remain in
the acceptance table. No production source, case controls or saved models
were modified for this plan.
