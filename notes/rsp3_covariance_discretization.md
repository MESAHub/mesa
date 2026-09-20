# RSP3 moment closure and discrete update

2026-09-20. Mathematical proposal and executable Python reference checks.
No production Fortran changes or new MESA runs are included in this stage.
The run evidence and published coefficient audit are in
[rsp3_covariance_closure.md](rsp3_covariance_closure.md).

## Scope and status

The target supports local and nonlocal convection, long evolutionary steps,
pulsations, and moving convective boundaries. Keep w, Pi and Phi on faces.
Keep the gas thermodynamic variables in cells and Lc at their interfaces.
No new moment variable or physical coefficient is needed for the candidate.

- [x] Check the published normalization and retain its distinction from this proposal.
- [x] Derive the local covariance identity exactly.
- [x] Derive a pressure update paired with the actual face energy work.
- [x] Verify residual derivatives independently in Python.
- [x] Test stable and unstable implicit local steps, including zero moments.
- [x] Establish the common transport and energy accounting requirements.
- [x] Test a boundary moving across several unequal face volumes in one step.
- [x] Test the large-step residual normalization with fixed reference scales.
- [x] Separate the physical closure, time discretization and nonlinear solver changes.
- [x] Map the proposal onto the existing MESA routines and validation gates.
- [ ] Translate the accepted equations, predictor and checks into MESA.
- [ ] Verify the native complete Jacobian, mesh paths, energy balance and LNA.
- [ ] Compare stellar evolution and pulsation with timestep and mesh refinement.

The user clarified that both local and nonlocal operation are required.
Enhanced stable-layer dissipation is out of scope for this implementation.
Do not add beta, a buoyancy dissipation length or a Brunt-frequency decay rate.

## Proposed implementation contract

The mathematical reference is sufficiently developed to start the native
implementation design. Further investigation should address the specific
integration checks below. The Python results do not establish stellar
convergence, calibration or pulsation growth rates.

Keep the three unknowns w, Pi and Phi on faces, with e_t=w^2. Their
covariance condition is part of the model: a finite convective flux cannot
coexist with zero entropy variance. Preserve the existing Y-based entropy
driving, EOS, mixing length and face luminosity relation. Do not replace
that driving with a difference of neighboring entropy values.

| Part | Proposed change | Consequence |
| --- | --- | --- |
| Local closure | Pair Pi buoyancy production and decay with e_t and Phi | Changes the physical moment dynamics; not a recovered literature factor |
| Compression and viscosity | Pair Pi compression with turbulent pressure work; pair Eq with the stress work in Uq | Preserve the moment condition and mechanical energy accounting together |
| Nonlocal closure | Use the existing alfat conductance for all three moments | Alfat=0 is local; nonzero alfat transports e_t, Pi and Phi |
| Time discretization | Use the new state for all moment transport and match Lt in gas energy | Changes finite-step Lt damping; retain the existing heat-flux time weight for Lr and Lc |
| Nonlinear solution | Use a coupled admissible predictor and fixed rate normalization | Does not add a physical source or change the residual's exact root |
| Remeshing | Apply the same positive conservative map to e_t, Pi and Phi | Preserves admissibility when donor states are admissible |

Use the existing RSP2_3equation_flag. The proposed normalization keeps
RSP2_alfa_pi and RSP2_alfa_phi at one by default and needs no new physical
control. RSP2_alfa_pi acquires the paired normalization defined below;
its documentation and allowed range must change with the source. Keep
ordinary RSP2 and TDC source closures outside this change.

With RSP2_3equation_flag false, retain the current one-equation RSP2
formulation, including its installed face w layout, algebraic convective
flux, turbulent-energy source, Lt time weighting, viscosity and predictor.
The paired moment decay, new moment transport, Lt-only gas-energy correction,
paired RSP3 viscous work and moment solver changes require the three-equation
flag. A shared routine edit must preserve the flag-off residuals and AD
partials to roundoff. Include local and nonzero-alfat flag-off comparisons.
No return to the older cell w implementation is planned.

The user chose to retain RSP2_alfa_pi at one. At that value its multiplier
is redundant with the paired decay coefficient. Values above one add
decorrelation; values below one can drive an admissible state outside the
covariance bound and are excluded from the proposed normalization.

Checkpoint the current face implementation before changing the closure.
Include its source, matching manuscripts and focused derivation/check scripts.
Exclude user test-case settings, generated run output and unrelated research.
Record the checkpoint as the current implementation with known moment-domain
failures, not as validation of the corrected closure. The first implementation
tasks are the predictor damping bound and the complete trial/acceptance path
specified below; they do not require a further general literature survey.

Do not prescribe Lc=0 in every stable layer. A signed transient flux is
allowed by this model. The objective is a consistent evolution of its
moments, not the elimination of every stable turbulent tail.

## Definitions and the physical change

Use the existing definitions

```math
e_t=w^2=\tfrac12\langle|\boldsymbol{u}'|^2\rangle,\qquad
\Pi=\langle u'_r s'\rangle,\qquad \Phi=\langle s'^2\rangle.
```

The full Phi here is twice the variance variable in Braun et al. (2026).
The current Fortran includes their isotropy and normalization factors.
The proposal below changes the reduced closure; it is not a transcription
repair to those published coefficients.

For fixed isotropy, admissible moments satisfy

```math
w\geq0,\quad \Phi\geq0,\quad
\Pi^2\leq\tfrac23w^2\Phi.
```

Both signs of Pi and countergradient heat flow remain possible. The source
uses the existing mapped entropy gradient from Y and the selected temperature
equation. The buoyancy coefficient is the existing code quantity

```math
\texttt{buoyancy}=-\frac{1}{\rho}\frac{\partial P}{\partial r}
                         \frac{\chi_T}{\chi_\rho c_p}.
```

Hold the energy source `buoyancy*Pi` and the Pi entropy source
`(2/3)*w^2*(-ds/dr)` fixed. Within an ansatz linear in the moments at fixed
gas state, a constant multiplier of the Pi buoyancy term must be one third
to cancel buoyancy production in the covariance determinant for both signs
of Pi. The exact identity before selecting this multiplier is recorded in
`symbolic_closure_checks.json`.

This can be interpreted as modelling the pressure/entropy correlation omitted
in the published reduced equations. It is a phenomenological closure choice.
The covariance proof does not establish its calibration, its superiority to
a model evolving the radial stress, or agreement with stellar growth rates.

Set

```math
C_D=\tfrac83\sqrt{\tfrac23}\,\texttt{RSP2_alfad},\qquad
C_\Phi=4\sqrt{\tfrac23}\,\texttt{RSP2_alfa_phi},\qquad
C_\Pi=\tfrac12(C_D+C_\Phi)\,\texttt{RSP2_alfa_pi}.
```

Require nonnegative variance decay rates and RSP2_alfa_pi >= 1 for this
interpretation. The default multipliers remain one. Unit controls preserve
the previous uniform stationary convective solution when radiation, work,
viscosity and transport vanish. Other coefficients and transient solutions
are not claimed to be unchanged.

The local reaction terms, before pressure work and transport, are

```math
\begin{aligned}
\dot e_t&=\texttt{buoyancy}\,\Pi-C_D\frac{w}{\Lambda}e_t+E_q,\\
\dot\Pi&=-\tfrac23e_t\frac{\partial s}{\partial r}
 +\tfrac13\texttt{buoyancy}\,\Phi
 -(C_\Pi w/\Lambda+\tau_{\rm rad}^{-1})\Pi,\\
\dot\Phi&=-2\Pi\frac{\partial s}{\partial r}
 -(C_\Phi w/\Lambda+2\tau_{\rm rad}^{-1})\Phi.
\end{aligned}
```

At the default Pi multiplier their covariance determinant obeys

```math
\frac{d}{dt}\left(\tfrac23e_t\Phi-\Pi^2\right)
=-\left[(C_D+C_\Phi)w/\Lambda+2/\tau_{\rm rad}\right]
  \left(\tfrac23e_t\Phi-\Pi^2\right)+\tfrac23E_q\Phi.
```

An additional nonnegative Pi decorrelation rate increases the right-hand
side by a nonnegative multiple of Pi squared. Nonnegative Eq is a condition
of this statement. The current new-strain times centered-strain work does
not guarantee it during velocity reversal. A dissipative revision must pair
the stress in Uq with the same work strain in Eq and preserve their discrete
mechanical-energy identity. It cannot clip Eq independently.

## Pressure work at finite timestep

`star_utils:calc_Ptrb_work_face` gives the face energy increment

```math
W_{P,f}=\sum_{j=k-1}^{k}\frac{\alpha_p\Delta m_j}{3\Delta m_f}
 [\theta_P\rho_j e_{t,f}+(1-\theta_P)\rho_{j,0}e_{t,f,0}]
 (\rho_j^{-1}-\rho_{j,0}^{-1}),\qquad
\Delta m_f=\tfrac12(\Delta m_{k-1}+\Delta m_k).
```

Write this exactly as `work_new*e_t + work_start*e_t_start`, where the
two work coefficients are AD quantities, not extra solver variables:

```math
\begin{aligned}
\texttt{work_new}&=\frac{\alpha_p\theta_P}{3\Delta m_f}
       \sum_j\Delta m_j(1-\rho_j/\rho_{j,0}),\\
\texttt{work_start}&=\frac{\alpha_p(1-\theta_P)}{3\Delta m_f}
       \sum_j\Delta m_j(\rho_{j,0}/\rho_j-1).
\end{aligned}
```

The energy residual keeps this exact work. A compatible Pi residual adds

```math
W_{\Pi,f}=\tfrac12\texttt{work_new}\,\Pi_f
 +\frac{\texttt{work_start}}
 {1+\sqrt{1-\texttt{work_start}}}\Pi_{f,0}.
```

The quotient evaluates `1-sqrt(1-work_start)` without subtracting nearby
numbers. It reduces the old covariance by the square root of the corresponding
old kinetic-variance factor. Merely adding half the new work while retaining
the full old Pi can violate the covariance bound.

This candidate requires `work_start < 1`. Removing more than the available
old kinetic energy through the explicitly weighted pressure work already
invalidates its energy-only update. A stellar step with such a density change
must be retried; arbitrary timesteps cannot override an inconsistent prescribed
hydrodynamic increment. The usual evolutionary small fractional density change
does not require resolving the convective turnover time.

In the continuum limit, W_Pi/dt is alpha_p*Pi*div(u)/3. Linearizing the actual
work about a stationary background gives the LNA inertia

```math
\delta\Pi_f-
 \frac{\alpha_p\Pi_f}{6\Delta m_f}
 (\Delta m_{k-1}\delta\ln\rho_{k-1}+\Delta m_k\delta\ln\rho_k).
```

This expression is independent of theta_P. Pi work belongs in the LNA mass
matrix, like the existing turbulent-energy pressure inertia. An unrelated
radial strain term must not be added a second time.

## Nonlocal transport and the heat flux

For nonzero alfat, use the same conductance for e_t, Pi and Phi, with the
current `compute_Lt_center` geometry:

```math
F_{a,k}=-\frac{\alpha_t(4\pi r_{c,k}^2\rho_k)^2\Lambda_{c,k}
                       (w_k+w_{k+1})/2}{\Delta m_k}
                  (a_k-a_{k+1}),\qquad a=e_t,\Pi,\Phi.
```

These are moment transport fluxes at gas cell centers. The mean enthalpy flux
remains `Lc=4*pi*r^2*rho_face*T_face*Pi` at gas interfaces. Use the same
boundary mask, positive face volumes and new-state conductance for all three
moments. Common implicit transport is a positive mixing of covariance matrices
and conserves their mass integrals with closed boundaries. Alfat=0 removes
all three transport terms and leaves the same local reaction model.

The existing energy transport is theta weighted by the luminosity control.
Unrestricted centered diffusion cannot guarantee nonnegative variance. In a
three-volume test with initial Phi=(0,1,0), unit conductances, theta=1/2 and
dt=10, it gives Phi=(0.625,-0.25,0.625). Fully implicit transport gives
(0.32258,0.35484,0.32258). All turbulent energies are positive in this test.

The preferred candidate for long evolutionary timesteps therefore makes
moment transport implicit. Preserve theta_L for radiative and convective
heat flow, and account for the different Lt weighting explicitly in the
total-energy row:

```math
L_{\rm used}=\theta_L L+(1-\theta_L)L_0
                 +(1-\theta_L)(L_t-L_{t,0}).
```

Here `L=Lr+Lc+Lt` remains the instantaneous total luminosity. This makes
`L_used=theta_L*(Lr+Lc)+(1-theta_L)*(Lr_0+Lc_0)+Lt`. It does not make
the convective heat flux fully implicit. Both appearances of turbulent
transport must use this same time level, or the thermal gas absorbs an
unintended difference. The extra inner Lt partial can reach k+2 and must
use the existing extended derivative storage.

This is a numerical change for Lt and can affect finite-step pulsation
damping. Its error must decrease with timestep refinement. It is preferable
to claiming that centered transport is positive for arbitrary steps.

## Residuals and the nonlinear solution

Before the existing fixed row normalization, the complete candidate rows are

```math
\begin{aligned}
R_{e,f}&=e_{t,f}-e_{t,f,0}+W_{P,f}
       -\Delta t\,\dot e_{t,f}
       +\Delta t(F_{e,k-1}-F_{e,k})/\Delta m_f,\\
R_{\Pi,f}&=\Pi_f-\Pi_{f,0}+W_{\Pi,f}
       -\Delta t\,\dot\Pi_f
       +\Delta t(F_{\Pi,k-1}-F_{\Pi,k})/\Delta m_f,\\
R_{\Phi,f}&=\Phi_f-\Phi_{f,0}
       -\Delta t\,\dot\Phi_f
       +\Delta t(F_{\Phi,k-1}-F_{\Phi,k})/\Delta m_f.
\end{aligned}
```

The dotted terms are only the local reactions written above. MESA retains w
as the unknown and its existing divided energy row where applicable.

Phi appears algebraically in these rows. A negative intermediate Newton Phi
need not make their evaluation undefined. It must not be accepted as a stellar
state. Removing trial clipping therefore requires changing both `Bdomain`
and the trial-variable domain check, plus an acceptance check covering all
three moments. `hydro_vars:unpack_xh` can retain checks on accepted/loaded
states. The line search must be audited for any call to that routine.

The admissible set is convex in (e_t,Pi,Phi), but not in (w,Pi,Phi). For
example, a straight midpoint between zero moments and an admissible fully
correlated state can be inadmissible in the latter coordinates. Enforcing
the full bound only at the endpoints of a w line search is insufficient.

Ordinary Newton and timestep continuation alone fail some large-step,
moving-boundary reference cases. A useful fallback is a coupled moment
predictor at fixed gas state. Freeze the current w-dependent coefficients
and solve the three moment rows together with a positive diagonal iteration
term added to both sides. The fixed point satisfies the original residual.
Use sufficient iteration damping for a positive covariance resolvent.
This preserves the predictor moments and avoids selecting an unphysical
large-step branch. It introduces no physical relaxation term or new timestep.

The reference implementation uses this predictor only where the first solve
fails; the stellar implementation should also reuse a good existing state.
Up to 200 small moment solves were required in one
12-face experiment. Do not claim that this cost or count equals MESA's global
Newton iterations. A band solve and suitable acceleration still need evaluation.

### Residual scale from short to long timesteps

The increment residual `Phi-Phi_start-dt*Phi_rhs` multiplies the floating-point
source cancellation error by dt. Normalizing it only by a fixed variance
scale therefore becomes unnecessarily restrictive in the stationary limit.
This is separate from the original clipping failure at Phi=0 and finite Pi.

Use the equivalent rate residual and a positive reference scale containing
both storage and source rates. For the variance row, one reference choice is

```math
\frac{(\Phi-\Phi_0)/\Delta t-\dot\Phi+
                     (F_{\Phi,k-1}-F_{\Phi,k})/\Delta m_f}
 {\Phi_{\rm scale}/\Delta t+
  2|\partial s/\partial r|_{\rm ref}\Pi_{\rm scale}+
  (C_\Phi w_{\rm ref}/\Lambda_{\rm ref}+2/\tau_{{\rm rad},{\rm ref}})\Phi_{\rm scale}
  +\text{transport rate reference}}.
```

The Pi row uses its corresponding production, decay and transport references.
Transport references sum the absolute conductance contributions multiplying
the local and neighboring moment scales. Freeze these references for a given
nonlinear solve. The Jacobian is then the physical residual derivative divided
by that fixed scale. Keep the variable correction scales in moment units;
they must not be replaced by these rate scales.

Short steps recover the state-increment normalization. Long steps measure
the relative balance of physical source and transport rates. No source is
removed and the exact solution is unchanged. Acceptance still requires the
physical moment constraint and the usual correction tests. This scaling
must not be used to accept an inconsistent moment state.

With this normalization, all 48 moving-boundary cases converge to a maximum
scaled rate residual of 3.93e-13, with nonnegative variances and covariance
ratios at most 1+2.3e-16. Fifteen cases use the positive predictor; the largest
predictor count is 200. The final generic nonlinear solve uses at most 77
function evaluations. These counts are reference-solver results, not predicted
MESA iteration counts.

## Verification record and remaining cases

The scripts and machine-readable results are in
`output/review/rsp3_work2_moments_20260920`.

| Check | Result and scope |
| --- | --- |
| Symbolic covariance and pressure identities | Exact SymPy identities, including Pi LNA inertia |
| Pressure rewrite | 1000 unequal-mass examples, maximum scaled error 1.32e-13 |
| Local nonlinear implicit steps | 800 stable/unstable cases; includes 40 zero-energy and 40 zero-variance starts; all admissible |
| Residual derivatives | 15000 derivatives, maximum scaled error 3.30e-15 against complex step; not native MESA partials |
| Energy pairing | Exact heat-flux time-weight identity; 1000 projection and 2000 viscous-work incidence checks |
| Common implicit remap/transport | Unequal mass tests preserve covariance and conserve integrated moments |
| Centered transport | Explicit negative-variance counterexample; 1000 tests verify the sufficient timestep restriction |
| Moving boundary | 48 nonlinear 12-face cases, alfat=0,0.1,1; boundary moves from face 4 to 1,4,9,12; dt=1e-3,1,1e3,1e6 |

The moving-boundary tests prescribe the mean stratification. They do not
solve stellar hydrostatic balance, luminosity, composition or mesh adaptation.
Generic Newton failed 13 cases; ordinary timestep continuation left six
failures. The coupled positive predictor recovered those six without changing
the physical equations, starting states or timesteps. All 48 final states
were admissible to roundoff. The largest increment-scaled residual was
7.28e-8. The subsequent rate normalization lowers the maximum scaled rate
residual to 3.93e-13 without changing the physical equations. The native
reference scales and Jacobian still require their own checks.

The remaining validation matrix must cover:

- Exact zero moments, finite Phi at zero w, either sign of Pi, both gradient
  signs, rapid sign reversal, and admissible versus invalid restart states.
- Zero and nonzero radiation, turbulent pressure, viscosity and transport.
- v_flag and u_flag, dedt and eps_grav, simple and conservative work,
  pressure reconstruction, pressure/heat time weights and moving boundaries.
- Every temperature-gradient form, dynamical gradL and composition gradients,
  with the same Y-based entropy driving and its full derivatives.
- Ordinary remeshing, split/merge AMR and envelope remeshing; unequal cells,
  forced non-turbulent boundaries, conservation and covariance after remap.
- Full native Jacobians including k+2 terms, both linear solver paths and
  the LNA mass and source matrices.
- Timestep/mesh refinement of relaxation rates, equilibrium structure,
  penetration depth, pulsation phase/growth and total-energy errors.

## Native implementation sequence

The entries below describe planned changes, not completed Fortran edits.
Keep the equations in the existing modules. Introduce a helper only where
it shares an actual coefficient or operator between consumers.

| Location | Planned change and required check |
| --- | --- |
| `hydro_rsp2:rsp2_moment_rhs` | Replace the Pi buoyancy and decay coefficients. Remove the independent radial strain when the paired pressure increment is added. Reuse the current AD entropy gradient, EOS and radiation rate. |
| `star_utils:calc_Ptrb_work_face` and `hydro_rsp2:do1_rsp2_moment_eqns` | Expose the two coefficients of the existing pressure increment if that avoids duplicating its geometry. Add W_Pi with complete density derivatives. Check both boundaries and work_start < 1. |
| `hydro_rsp2` viscosity routines | Pair the RSP3 stress and work strain in Uq and Eq. Verify the native u/v mechanical-energy sums and both gas-energy work forms. A positive-square reference test alone does not verify the native geometry. |
| `hydro_rsp2:compute_Lt_center` and `rsp2_dLt_dm_face` | Share the existing conductance and closed-boundary mask with Pi and Phi transport. Include the derivatives of the conductance as well as the moment differences. Use the new state for RSP3 moment transport. |
| `hydro_energy:setup_dL_dm` | Add the Lt-only time-weight correction, using the actual local L_theta. Preserve the unshifted inner Lt derivative before shift_p1 and retain its k+2 contribution. Check dedt and eps_grav. |
| `hydro_rsp2:set_etrb_start_vars` and `do1_rsp2_moment_eqns` | Establish fixed rate references for each nonlinear solve. Keep Pi_scale and Phi_scale in state units for correction tests. Define reference initialization and invalidation at startup, retry and remesh; avoid stale caches. |
| `hydro_rsp2:RSP2_adjust_vars_before_call_solver` | Replace the independent moment seed with a coupled predictor where needed. Solve for energy internally and return nonnegative w. Reuse a native band solver; do not copy the reference's dense eigenvalue calculation. |
| `solver_support:Bdomain`, `set_vars_for_solver`, `star_solver` and `hydro_vars:unpack_xh` | Trace every trial and acceptance path before changing clipping. Check the full moment condition on accepted states, including the last allowed iteration. Keep w trials in their domain. A straight line in w does not inherit convexity in energy. |
| `hydro_rsp2:init_rsp2_moments` and the existing common remap | Check exact zero states, admissible active seeds and identical moment weights in ordinary remesh, split/merge and envelope remesh. Do not silently relabel an invalid photo as an admissible initial state. |
| `star_LNA_support:assemble_rsp2_moment_rows` and `star_LNA_turbulence_closures` | Reuse the changed physical sources and moment transport in A. Add the Pi density inertia to B. Predictor damping and residual reference scales are numerical devices and must not become physical LNA terms. |
| `controls_dev.defaults` and `set_flags` | Document the changed alfa_pi normalization, constrain the allowed decay rates and document all three alfat fluxes. Retain the unit default multipliers. |

Before coding the predictor, specify a sufficient inexpensive damping bound
for the native band system, its stopping test and its failure return. The
reference's successful 12-face solves are not a cost or convergence guarantee
for a full stellar mesh. Preserve the original start state and timestep
throughout predictor iteration. Do not accept a predictor as the stellar
solution without checking the complete coupled residual.

Before changing domain handling, trace calls from the line search to
unpack_xh and identify a final acceptance check that cannot be bypassed.
An algebraically evaluable negative Phi trial is different from accepting
negative variance. Test the full w, Pi, Phi path, including a zero-energy
start where the derivative of w^2 vanishes.

Implement and review the local source/work changes first, then transport
and its gas-energy accounting, then the predictor and solver integration.
Update LNA and initialization in the same implementation series before
stellar testing. Finish with remesh and control-off regression checks.
The existing manuscripts describe the installed equations; update them
after the new equations are in the source so that neither document claims
an unimplemented prescription is current.

Native validation requires source checks, complete AD partial checks and,
when authorized, builds and stellar runs. Assess moment admissibility,
residual convergence, total energy and pulsation growth separately. Enhanced
dissipation and its controls remain deferred.
