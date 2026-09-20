# RSP3 face profiles and local modes

Date: 2026-09-20. Saved-output analysis only. No additional MESA run,
compilation, case edit or physical closure change was made for this audit.

The face layout does not by itself remove the stable-layer moment pattern.
The earlier claim that placement explained the whole pattern was too broad.
The inner bumps and the selected extra LNA modes also occupy different layers.

## Data and reproducibility

The stopped isolated run is
`output/review/rsp3_face_w_run_20260920/alloc_fix`.
It has every accepted profile through model 1655 and the LNA at model 200.
The user's saved photos at models 1000, 2000 and 3000 were copied to
`output/review/rsp3_face_w_run_20260920/convection_audit` before analysis.
The script there, `analyze.py`, reads these files without invoking MESA.

The photo reader checks the record markers, version 22, nine hydro variables,
RSP3 flag and variable indices. At model 1000 it reproduces the independently
written profile density, temperature and radius within 2.3e-16 relative;
w and Y agree exactly. Pi agrees with w*PII_face within 1.3e-16 relative
to its maximum, and the reconstructed Lc agrees within 8.2e-16.
The user's model 1000 w is identical to the isolated run's model 1000 w.

For this case, with mass interpolation and no face reconstruction,

```math
L_{c,k}=4\pi r_k^2\rho_{\mathrm{face},k}T_{\mathrm{face},k}\Pi_k.
```

This is `hydro_rsp2:compute_Lc_terms`, using
`hydro_gradient_support:get_rsp2_face_eos`. The raw profile Lc is in erg/s.
The plots explicitly convert it to solar luminosities where indicated.
The inner core is excised; full radius means the whole modeled envelope.

Figures in the same directory:

- `moments_full_radius.png` and `.pdf`: Pi, Phi, Lc/Lsun and Lc/L over the envelope.
- `inner_pattern.png` and `.pdf`: the specific inner bumps in the user's screenshot.
- `inner_startup.png` and `.pdf`: their first appearance and subsequent remesh.
- `convection_model1000.png` and `.pdf`: w, moments, flux and saved D_mix.
- `convection_summary.json` and `inner_startup_budget.json`: numerical results.

## Inner pattern

The region is approximately 0.15 < r/R < 0.34. It is stable throughout
these snapshots. The strongest inner bump follows this sequence:

| Model | Zones | Maximum w in the region, cm/s | Event |
| --- | --- | --- | --- |
| 1 | 149 | 1.09e-7 | Initial accepted model |
| 2 | 350 | 6.58 | First ordinary mesh adjustment and subsequent hydro step |
| 100 | 477 | 6.39 | Before the envelope remesh |
| 101 | 350 | 2.85 | After the envelope remesh |
| 1000 | 350 | 2.84075 | After the LNA kick |
| 2000 | 350 | 2.84074 | Nearly unchanged |
| 3000 | 350 | 2.84083 | Nearly unchanged |

Thus the bumps predate both the model-200 LNA and its velocity kick.
Y and radial velocity also acquire structure during the first mesh change
and hydro step. Accepted profiles do not separate the immediate remap from
the subsequent Newton solve; they do not prove that remapping directly
deposited turbulent energy. At the model-2 peak, saved w_start is zero.

The energy equation identifies the actual source of the increase. In this
case alfat=alfap=alfar=0 and v_flag is active, so the face energy balance is

```math
\frac{w_k^2-w_{k,\mathrm{start}}^2}{\Delta t}
=\mathrm{SOURCE}_k-\mathrm{DAMP}_k+E_{q,k}^{\mathrm{face}}.
```

For v_flag the saved Eq is a cell value. The script first divides it by
(w_k+w_{k+1})/2, then applies the mass weights and multiplies by w_k,
matching `hydro_rsp2:compute_Eq_div_w_face`. It does not compare a cell Eq
directly to the face equation.

At model 2, face 326, all rates in erg/(g s):

| Quantity | Value |
| --- | --- |
| Face Eq | +5.77524374e-4 |
| SOURCE | -1.44356388e-4 |
| DAMP | 4.04892379e-8 |
| Sum of RHS | 4.33127497e-4 |
| (w^2-w_start^2)/dt | 4.33127433e-4 |

The increase is supplied by eddy-viscous heating during the startup hydro
response, while buoyancy opposes it. This identifies the energy path, not
a verdict that the startup shear or remeshed equilibrium is acceptable.
The zero-state predictor also explicitly permits a positive w guess from
Eq/w in `RSP2_adjust_vars_before_call_solver`.

At model 1000 the peak's w^2/DAMP is 81.47 years. Phi's direct turnover
decay is comparably slow; radiative cooling is disabled. Therefore tiny
moments can retain startup structure for many pulsation cycles.
The implemented local Pi equation also permits approximate cancellation
between its negative entropy-gradient production and positive buoyancy
production from Phi. Using the profile Brunt frequency and the discrete
buoyancy coefficient, the peak's Phi is within 0.3 percent of that balance.
There is no Pi/Phi transport to smooth differences between adjacent faces.

In this inner region, w and Phi barely change between models 1000 and 3000.
Pi and Lc change sign, but max |Lc/L| is respectively 1.46e-11, 1.05e-11
and 9.77e-12. The visible w/Phi structure is therefore a persistent spatial
imprint; these snapshots do not show large oscillatory heat transport there.
It still warrants investigation and can keep moment rows active in LNA.

`mix_info:set_mixing_info` computes D_mix = sqrt(2/3)*w*Lambda/3 and labels
stable nonzero turbulence as overshoot. This label does not establish
nonlocal transport. Saved D_mix is from the mixing update and need not
equal a fresh evaluation using the final Newton state on a rapidly changing
step; the plots use the actual saved diagnostic.

## Extra LNA modes

Fourteen of the fifteen selected periods match a local moment oscillator
estimate to within 0.7 percent, evaluated at the maximum |delta w| face.
Those maxima are at faces 176--185 and 199--202, around r/R=0.91--0.97,
not at the deep inner bumps discussed above.

For fixed gas coefficients, omitting damping, strain and transport, the
implemented local moment equations are

```math
\dot e_t=-\frac{1}{\rho}\frac{dP}{dr}
 \frac{\chi_T}{\chi_\rho c_p}\Pi,\qquad e_t=w^2,
```

```math
\dot\Pi=-\frac23e_t\frac{ds}{dr}
 -\frac{1}{\rho}\frac{dP}{dr}\frac{\chi_T}{\chi_\rho c_p}\Phi,
\qquad \dot\Phi=-2\Pi\frac{ds}{dr}.
```

They give approximately

```math
\ddot\Pi+\frac83N^2\Pi=0,\qquad
P_{\mathrm{local}}\simeq\frac{2\pi}{\sqrt{8/3}\,N}.
```

The comparison uses saved Brunt N2, rather than claiming exact equality to
the discrete gradient mapping. Code: `rsp2_moment_rhs`, `rsp2_moment_source`
and their LNA assembly. Exact dormant zero moments have algebraic rows;
small nonzero moments retain the dynamic rows. Small amplitude therefore
does not imply absence of a local eigenfrequency.

This is evidence for a local moment branch of the implemented closure,
not fourteen additional acoustic radial modes or proof that each computed
root is an accurate, physical stellar mode. Mode 13, P=0.70188259 days,
is the acoustic fundamental candidate: its broad displacement and residual
1.14e-10 distinguish it from this sequence. See the parent directory's
`lna_mode_audit.png` and `lna_mode_shape_audit.json`.

A separate accuracy issue remains: `select_modes` uses
min(user_max_residual, 0.1) as the refinement target. With the case's
max residual of 1, reported convergence can mean only residual <= 0.1.
The componentwise residual is itself bounded by 1. These settings cannot
certify the less accurate roots in the list. No tolerance was changed here.

## Remaining closure and derivative checks

The saved stable faces satisfy Pi^2 <= (2/3)*w^2*Phi in all three photos.
Small violations occur in unstable faces: the largest ratio is 1.0228 at
model 1000. This is separate from the inner stable pattern and does not
establish its cause. The previously derived covariance correction remains
unimplemented; see `rsp3_covariance_closure.md`.

The k+2 derivatives already have storage and an expanded band in
`star_solver:solve_with_banded_solver` (d_hydro_d_p2; band_ku=3*nvar when
needed). They arise from inner-face viscous work in conservative v_flag
energy, inner-face turbulent pressure work, and u_flag Riemann pressure
and work using the adjacent cell turbulent pressure. The local moment
sources and Lt transport remain adjacent-zone stencils. The present case's
v_flag, simple work and alfap=0 do not activate these extra terms.

Next discriminating work: inspect the immediate ordinary-remesh thermal
and pressure balance and the ensuing Eq production, separately from the
long stable-layer decay time. A controlled run is required to verify a
remesh change. Neither clipping the tails nor modifying LNA selection
would fix their creation. A buoyancy decay term is a physical closure
change and should not be presented as a demonstrated numerical repair.

The user suggested relaxation into RSP3 or larger initial timesteps. This
is consistent with a startup imprint, but has not been tested here. Settle
the background on the final mesh before assessing the moments. Larger
implicit steps can damp oscillatory moment components numerically; they
do not prove that a nearly stationary w/Phi component has physically
decayed. Restore the intended pulsation timesteps before measuring growth
or deciding whether a closure or LNA change is required.
