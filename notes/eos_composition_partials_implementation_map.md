---
title: "EOS Composition Partials: Implementation Map"
subtitle: "Skye, eosDT, Brunt, and implicit diffusion follow-up"
date: "2026-05-14"
geometry: margin=0.55in
fontsize: 10pt
colorlinks: true
linkcolor: blue
urlcolor: blue
header-includes:
  - \usepackage[most]{tcolorbox}
  - \usepackage{xcolor}
  - \usepackage{booktabs}
  - \usepackage{array}
  - \usepackage{enumitem}
  - \usepackage{amsmath}
  - \usepackage{amssymb}
  - \usepackage{fvextra}
  - \DefineVerbatimEnvironment{Highlighting}{Verbatim}{breaklines,commandchars=\\\{\}}
  - \definecolor{currentbg}{HTML}{EAF1FB}
  - \definecolor{currentframe}{HTML}{2F5E9E}
  - \definecolor{newbg}{HTML}{EAF7EA}
  - \definecolor{newframe}{HTML}{2B7A3D}
  - \definecolor{starbg}{HTML}{FFF1DC}
  - \definecolor{starframe}{HTML}{B86B00}
  - \definecolor{verifybg}{HTML}{F3EAFE}
  - \definecolor{verifyframe}{HTML}{6F42C1}
  - \definecolor{riskbg}{HTML}{FDECEC}
  - \definecolor{riskframe}{HTML}{B3261E}
  - \definecolor{decisionbg}{HTML}{FFF9D7}
  - \definecolor{decisionframe}{HTML}{9A7A00}
  - \newtcolorbox{currentbox}{colback=currentbg,colframe=currentframe,boxrule=0.7pt,arc=1mm,left=2mm,right=2mm,top=1mm,bottom=1mm}
  - \newtcolorbox{newimplbox}{colback=newbg,colframe=newframe,boxrule=0.7pt,arc=1mm,left=2mm,right=2mm,top=1mm,bottom=1mm}
  - \newtcolorbox{starbox}{colback=starbg,colframe=starframe,boxrule=0.7pt,arc=1mm,left=2mm,right=2mm,top=1mm,bottom=1mm}
  - \newtcolorbox{verifybox}{colback=verifybg,colframe=verifyframe,boxrule=0.7pt,arc=1mm,left=2mm,right=2mm,top=1mm,bottom=1mm}
  - \newtcolorbox{riskbox}{colback=riskbg,colframe=riskframe,boxrule=0.7pt,arc=1mm,left=2mm,right=2mm,top=1mm,bottom=1mm}
  - \newtcolorbox{decisionbox}{colback=decisionbg,colframe=decisionframe,boxrule=0.7pt,arc=1mm,left=2mm,right=2mm,top=1mm,bottom=1mm}
---

This document is a rendered companion to
`notes/eos_composition_partials_plan.md`.  It separates what MESA already
does from what must be added, and it flags implementation effort.  No MESA
compile, MESA run, test suite, or EOS table regeneration has been done.

# Legend

\begin{currentbox}
\textbf{Current MESA behavior.}  Existing code or equations already present
in the repository.  These are source-level facts to preserve or use.
\end{currentbox}

\begin{newimplbox}
\textbf{New EOS implementation.}  Analytic EOS-side derivatives to add,
mostly in \texttt{eos/private}.  This is the main work needed to remove
finite-difference composition partials for Skye/eosDT-covered cells.
\end{newimplbox}

\begin{starbox}
\textbf{New or changed star-side use.}  Changes in \texttt{star/private}
that consume the new EOS partials, especially Brunt and possible coupled
implicit diffusion work.
\end{starbox}

\begin{verifybox}
\textbf{Verification only.}  Finite differences are acceptable here as
tests of analytic derivatives, but not as production replacement partials.
\end{verifybox}

\begin{decisionbox}
\textbf{Design decision.}  A choice that should be made explicitly before
implementation broadens.
\end{decisionbox}

\begin{riskbox}
\textbf{Risk or nonlinear boundary.}  Active-set changes, blend kinks,
phase branches, and API expansion points.
\end{riskbox}

# Effort Summary

| Area | Color | Effort | Why |
|---|---:|---:|---|
| Composition helper and moment derivatives | New EOS | S | Mostly wraps formulas already in `chem_lib`. |
| Skye ideal-ion partials | New EOS | M | Direct analytic chain rule through active-set number fractions. |
| Skye electron partials | New EOS | M | HELM derivative tables already contain the needed density derivatives. |
| Skye Coulomb partials | New EOS | L/XL | Many scalar composition inputs and phase/branch logic. |
| eosDT blend `dalpha/dX_i` | New EOS | M | Mechanical plumbing plus blend alpha derivatives. |
| Brunt chain-rule form | Star-side | M | Branch exists, but needs total-pressure conversion and analytic EOS partials. |
| Coupled implicit diffusion | Star-side | XL | Separate solver-structure project; EOS partials are prerequisite. |
| Rotation/other transport compatibility | Star-side | M | Default-off controls should preserve current operator-split workflows. |
| Finite-difference checks | Verification | M | Needed to validate; not production behavior. |

# Current Contract

\begin{currentbox}
\textbf{EOS public export.}

\[
\texttt{num\_eos\_d\_dxa\_results}=2
\]

Only the first two composition derivative rows are reliably returned to
\texttt{star}:

\[
\texttt{d\_dxa(i\_lnPgas,j)}
  =
  \left.\frac{\partial \ln P_{\rm gas}}{\partial X_j}\right|_{\rho,T},
\qquad
\texttt{d\_dxa(i\_lnE,j)}
  =
  \left.\frac{\partial \ln E}{\partial X_j}\right|_{\rho,T}.
\]

Internally \texttt{eos/public/eos\_lib.f90} already allocates
\texttt{num\_eos\_basic\_results} rows and then copies only the first two
rows back to \texttt{star}.
\end{currentbox}

\begin{currentbox}
\textbf{Composition normalization and gauge.}

MESA composition variables satisfy

\[
S_X \equiv \sum_i X_i = 1,
\qquad
\sum_i \delta X_i = 0
\]

for physical perturbations.  A species-basis derivative vector has gauge
freedom: adding a constant to all species partials does not change any
constrained derivative or any contraction with a zero-sum gradient.

\[
\sum_i
\left(
  \frac{\partial Q}{\partial X_i}+C
\right)
\delta X_i
=
\sum_i
\frac{\partial Q}{\partial X_i}\delta X_i.
\]

\texttt{star/private/star\_solver.f90} currently uses a sink projection:

\[
D_j^{\rm sink}Q
=
\frac{\partial Q}{\partial X_j}
-
\frac{\partial Q}{\partial X_s}.
\]
\end{currentbox}

# Composition Moments

\begin{currentbox}
\textbf{Existing equations in \texttt{chem/public/chem\_lib.f90}.}

Let

\[
Y_i=\frac{X_i}{A_i},\qquad
Y=\sum_i Y_i,\qquad
\bar A = \frac{S_X}{Y}.
\]

For a charge moment \(g_i\),

\[
M_g = \bar A\sum_i\frac{X_i}{A_i}g_i.
\]

MESA uses \(g_i=Z_i\), \(Z_i^2\), and \(\texttt{chem\_isos\%Z53}\) for
\(\bar Z\), \(\overline{Z^2}\), and \(\overline{Z^{5/3}}\).  The
implemented derivatives are

\[
\frac{\partial \bar A}{\partial X_j}
=
\frac{\bar A}{A_jS_X}(A_j-\bar A),
\]

\[
\frac{\partial M_g}{\partial X_j}
=
\frac{\bar A}{A_jS_X}(g_j-M_g),
\]

\[
m_c=\sum_i X_i\frac{W_i}{A_i},
\qquad
\frac{\partial m_c}{\partial X_j}
=
\frac{W_j}{A_j}-m_c.
\]
\end{currentbox}

\begin{newimplbox}
\textbf{Add reusable EOS helper.}

Add a private helper, probably \texttt{eos/private/eos\_composition\_partials.f90},
that returns

\[
\bar A_j,\quad \bar Z_j,\quad Y_{e,j},\quad
\overline{Z^2}_j,\quad \overline{Z^{5/3}}_j,\quad m_{c,j},
\]

plus active-set normalized number-fraction derivatives for Skye.

For complete ionization,

\[
Y_e=\frac{\bar Z}{\bar A}
    =\sum_i X_i\frac{Z_i}{A_i},
\qquad
Y_{e,j}=\frac{Z_j}{A_j}-Y_e.
\]

Effort: \textbf{S}.  These are direct formulas already consistent with
\texttt{chem\_lib}.
\end{newimplbox}

# Skye Analytic Partials

\begin{currentbox}
\textbf{Current Skye state.}

\texttt{eos/private/skye.f90:Get\_Skye\_EOS\_Results} calls
\texttt{skye\_eos(..., d\_dxa, ierr)} and then sets

\[
\texttt{d\_dxa}=0.
\]

The value path already uses \texttt{auto\_diff\_real\_2var\_order3} for
\(T\) and \(\rho\), but not for arbitrary composition variables.
\end{currentbox}

\begin{newimplbox}
\textbf{Core Skye strategy.}

Do not introduce global arbitrary-size autodiff.  For each species \(j\),
compute a companion object

\[
F_j
=
\left.\frac{\partial F}{\partial X_j}\right|_{\rho,T},
\]

where \(F_j\) itself is still an \texttt{auto\_diff\_real\_2var\_order3}
object in \(T\) and \(\rho\).

Skye's existing thermodynamics are

\[
S=-\frac{\partial F}{\partial T},\qquad
P_{\rm gas}=\rho^2\frac{\partial F}{\partial \rho},\qquad
E=F+TS.
\]

Therefore

\[
S_j=-\frac{\partial F_j}{\partial T},\qquad
P_{{\rm gas},j}=\rho^2\frac{\partial F_j}{\partial \rho},\qquad
E_j=F_j+TS_j.
\]

The packed outputs are

\[
\frac{\partial\ln S}{\partial X_j}=\frac{S_j}{S},\qquad
\frac{\partial\ln E}{\partial X_j}=\frac{E_j}{E},\qquad
\frac{\partial\ln P_{\rm gas}}{\partial X_j}
=\frac{P_{{\rm gas},j}}{P_{\rm gas}}.
\]

Effort: \textbf{M} for the framework, before Coulomb.
\end{newimplbox}

\begin{newimplbox}
\textbf{Ideal-ion term.}

Skye's ideal-ion free energy is

\[
F_{\rm ion}
=
\frac{k_BT}{m_u\bar A}\Phi,
\]

\[
\Phi
=
\sum_i y_i
\left[
\ln\left(\frac{y_i n}{n_{Q,i}}\right)-1
\right],
\qquad
n=\frac{\rho}{m_u\bar A},
\qquad
n_{Q,i}=n_QW_i^{3/2}.
\]

With active-set normalized number fractions \(y_i\),

\[
\Phi_j
=
\sum_i y_{i,j}
\ln\left(\frac{y_i n}{n_{Q,i}}\right)
-\frac{\bar A_j}{\bar A},
\]

\[
F_{{\rm ion},j}
=
\frac{k_BT}{m_u}
\left(
\frac{\Phi_j}{\bar A}
-
\frac{\Phi\bar A_j}{\bar A^2}
\right).
\]

Ion offsets use full-network number fractions:

\[
F_{\rm off}
=
\frac{\rm eV}{m_u}\sum_i \tilde y_i I(Z_i),
\qquad
F_{{\rm off},j}
=
\frac{\rm eV}{m_u}\sum_i \tilde y_{i,j} I(Z_i).
\]

Effort: \textbf{M}.  Main risk is keeping active-set and full-network
normalizations distinct.
\end{newimplbox}

\begin{newimplbox}
\textbf{Electron term.}

\texttt{compute\_ideal\_ele} evaluates the HELM free energy as

\[
D=Y_e\rho,\qquad
F_e=Y_e f(T,D).
\]

For

\[
f_{a,b}
=
\frac{\partial^{a+b}f}{\partial T^a\partial D^b},
\]

the stored coefficient is

\[
\frac{\partial^{a+b}F_e}{\partial T^a\partial\rho^b}
=
Y_e^{b+1}f_{a,b}.
\]

The composition derivative is

\[
\frac{\partial}{\partial X_j}
\left[
\frac{\partial^{a+b}F_e}{\partial T^a\partial\rho^b}
\right]
=
Y_{e,j}
\left[
(b+1)Y_e^b f_{a,b}
+\rho Y_e^{b+1} f_{a,b+1}
\right].
\]

Current Skye overwrites HELM \(\texttt{xnefer}\) with complete-ionization
\(\texttt{compute\_xne}\), so

\[
\ln{\rm free}_e
=
\ln Y_e,
\qquad
\frac{\partial\ln{\rm free}_e}{\partial X_j}
=
\frac{Y_{e,j}}{Y_e},
\]

away from the existing \(\max(10^{-99},x_{\rm nefer})\) clamp.

Effort: \textbf{M}.  Most table derivatives already exist.
\end{newimplbox}

\begin{newimplbox}
\textbf{Coulomb term.}

Treat Coulomb as a chain-rule problem over composition-dependent scalar
inputs:

\[
F_C=G(T,\rho,u_1,\ldots,u_m),
\qquad
F_{C,j}
=
\sum_a
\frac{\partial G}{\partial u_a}
\frac{\partial u_a}{\partial X_j}.
\]

The implementation should use a Skye-local companion type carrying
\((u,u_1,\ldots,u_N)\), where each entry is an
\texttt{auto\_diff\_real\_2var\_order3} in \(T,\rho\).  This avoids
hand-expanding every Coulomb fit.

Effort: \textbf{L/XL}.  This is the largest EOS-side piece because of
liquid/solid fits, mixing entropy, phase softmin, and branch behavior.
\end{newimplbox}

\begin{riskbox}
\textbf{Skye branch risks.}

Derivatives are only smooth while the active set is unchanged:

\[
X_i > \texttt{mass\_fraction\_limit\_for\_Skye}.
\]

At exact active-set boundaries, \(\min/\max\) boundaries, and phase branch
boundaries, use the derivative of the selected value branch and avoid
those points in finite-difference tests.
\end{riskbox}

# eosDT and Blends

\begin{currentbox}
\textbf{Current table mapping.}

\texttt{eos/private/eosdt\_eval.f90:get1\_for\_eosdt} maps table
derivatives by charge:

\[
\frac{\partial R}{\partial X_j}
=
\frac{\partial R}{\partial X_H}H_j
+
\frac{\partial R}{\partial Z}M_j,
\]

with

\[
H_j=
\begin{cases}
1,&Z_j=1,\\
0,&\hbox{otherwise},
\end{cases}
\qquad
M_j=
\begin{cases}
1,&Z_j>2,\\
0,&\hbox{otherwise}.
\end{cases}
\]

Helium gets zero table derivative because the tables use
\(Y=1-X_H-Z\).
\end{currentbox}

\begin{newimplbox}
\textbf{Blend derivative to add.}

For

\[
R=\alpha R_1+(1-\alpha)R_2,
\]

the complete composition derivative is

\[
R_j
=
\alpha R_{1,j}
+
(1-\alpha)R_{2,j}
+
\alpha_j(R_1-R_2).
\]

Current \texttt{Do\_Blend} has the first two terms only.  Add
\(\alpha_j\), including the existing quintic smoothing:

\[
\alpha=h(a_0)=10a_0^3-15a_0^4+6a_0^5,
\]

\[
\alpha_j=h'(a_0)a_{0,j},
\qquad
h'(a_0)=30(a_0-1)^2a_0^2.
\]

Effort: \textbf{M}.  The API plumbing is mechanical; alpha derivatives
need case-by-case audit.
\end{newimplbox}

# Star-Side Consumers

\begin{currentbox}
\textbf{Current Ledoux/MHM approximation.}

Current MESA computes the Brunt composition term in
\texttt{star/private/brunt.f90:do\_brunt\_B\_MHM\_form} by calling the EOS
twice at the same face \((\rho_f,T_f)\), but with neighboring cell
compositions:

\[
\ln P_1=\ln P_{\rm eos}(\rho_f,T_f,X_k),
\qquad
\ln P_2=\ln P_{\rm eos}(\rho_f,T_f,X_{k-1}).
\]

Then

\[
B_{\rm MHM}
=
\frac{\ln P_1-\ln P_2}
     {(\ln P_{{\rm eos},k-1}-\ln P_{{\rm eos},k})\chi_{T,f}}.
\]

Linearizing the numerator gives

\[
\ln P_1-\ln P_2
\approx
\sum_j\chi_{X_j,f}(X_{j,k}-X_{j,k-1}),
\]

so

\[
B_{\rm MHM}
\approx
-\frac{1}{\chi_{T,f}}
\sum_j\chi_{X_j,f}
\frac{X_{j,k-1}-X_{j,k}}
     {\ln P_{{\rm eos},k-1}-\ln P_{{\rm eos},k}}.
\]

That is the same first-order Ledoux chain-rule structure, but it gets the
EOS composition contraction from a finite pressure difference.  After
\texttt{brunt\_B} is computed, \texttt{hydro\_vars.f90} stores it in
\texttt{gradL\_composition\_term}, and \texttt{turb\_support.f90} uses

\[
\nabla_L=\nabla_{\rm ad}+\texttt{gradL\_composition\_term}.
\]
\end{currentbox}

\begin{starbox}
\textbf{Brunt chain-rule form.}

The branch \texttt{origin/Brunt\_B\_from\_eos\_partials} wants

\[
B_{\rm comp}
=
-\frac{1}{\chi_T}
\sum_j \chi_{X_j}\frac{dX_j}{d\ln P}.
\]

The EOS partial itself must be returned by the EOS at the face state:

\[
g_{j,f}
=
\left.
\frac{\partial\ln P_{\rm gas}}{\partial X_j}
\right|_{\rho_f,T_f,X_f}.
\]

This is the part we are trying to fix.  It should not be approximated by
finite-differencing EOS calls with neighboring cell compositions.  Since
MESA currently returns gas-pressure partials, the total-pressure
conversion must then be applied at the same face:

\[
\chi_{X_j}
=
\frac{P_{\rm gas}}{P_{\rm eos}}
g_{j,f},
\]

\[
\chi_T
=
\frac{P_{\rm gas}}{P_{\rm eos}}
\left.
\frac{\partial\ln P_{\rm gas}}{\partial\ln T}
\right|_{\rho,X}
+
\frac{4P_{\rm rad}}{P_{\rm eos}}.
\]

The remaining finite difference is the spatial composition gradient, not
an EOS composition partial:

\[
\frac{dX_j}{d\ln P}
\approx
\frac{X_{j,k-1}-X_{j,k}}
     {\ln P_{{\rm eos},k-1}-\ln P_{{\rm eos},k}}.
\]

Effort: \textbf{M}, after analytic EOS partials exist.
\end{starbox}

\begin{starbox}
\textbf{Implicit diffusion coupling.}

Current element diffusion is implicit inside
\texttt{star/private/diffusion.f90:do\_solve\_diffusion}, but it is
operator-split from the hydro Newton solve.  MESA already has one hydro
composition equation per isotope in
\texttt{star/private/hydro\_chem\_eqns.f90:do1\_chem\_eqns}; coupled
microscopic diffusion should add a rate to that equation, not create a
second species equation.

In current term names,

\[
\mathrm{dxdt\_actual}_{j,k}
=
\frac{X^{n+1}_{j,k}-X^n_{j,k}}{\Delta t},
\qquad
\mathrm{dxdt\_expected}_{j,k}
=
\mathrm{dxdt\_mix}_{j,k}
+
\mathrm{dxdt\_nuc}_{j,k},
\]

\[
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
\]

The coupled diffusion term should be

\[
\mathrm{dxdt\_expected}_{j,k}
=
\mathrm{dxdt\_mix}_{j,k}
+
\mathrm{dxdt\_nuc}_{j,k}
+
\mathrm{dxdt\_diff}_{j,k},
\]

\[
\mathrm{dxdt\_diff}_{j,k}
=
-
\frac{\Phi_{j,k+1/2}-\Phi_{j,k-1/2}}{\Delta m_k}.
\]

To allow some transport to remain split while another part is coupled,
accumulate the coefficient in two pieces:

\[
D_{\rm mix}
=
D_{\rm mix}^{\rm explicit}
+
D_{\rm mix}^{\rm implicit}.
\]

Here explicit means outside the main hydro Newton block, or included as a
frozen source with no new Jacobian.  The implicit part supplies the new
rate and Jacobian:

\[
\mathrm{dxdt\_diff,implicit}_{j,k}
=
-
\frac{
\Phi^{\rm implicit}_{j,k+1/2}
-
\Phi^{\rm implicit}_{j,k-1/2}
}{\Delta m_k},
\qquad
\frac{\partial R^X_{j,k}}{\partial X_{\ell,m}}
\supset
\frac{1}{\mathrm{eqn\_scale}_{j,k}}
\frac{\partial
\mathrm{dxdt\_diff,implicit}_{j,k}}
{\partial X_{\ell,m}}.
\]

The mass-fraction constraint can be handled by a sink species or by

\[
R^\Sigma_k=\sum_jX_{j,k}-1=0.
\]

Effort: \textbf{XL}.  EOS composition partials are prerequisite, not the
whole implicit-diffusion implementation.
\end{starbox}

\begin{starbox}
\textbf{Rotation and other transport should stay opt-in.}

Default behavior should not change: \texttt{do\_element\_diffusion} keeps
using the current operator-split diffusion solver, and rotation/ordinary
mixing keep entering through \texttt{s\%dxdt\_mix}.  Rotation modifies
\texttt{s\%D\_mix} through
\texttt{star/private/mix\_info.f90:update\_rotation\_mixing\_info}:

\[
D_{\rm mix}=D_{\rm mix,nonrot}+D_{\rm mix,rot}.
\]

A low-risk control structure is:

\begin{tabular}{ll}
\texttt{do\_element\_diffusion} & current operator-split path\\
\texttt{use\_coupled\_element\_diffusion} & new default-off Newton path\\
\texttt{use\_coupled\_diffusion\_rotation\_terms} & default off\\
\texttt{use\_coupled\_diffusion\_other\_D\_mix} & default off
\end{tabular}

Internally:

\[
\begin{aligned}
D_{\rm mix,explicit} &=
D_{\rm mix,nonrot}+D_{\rm mix,rot}+\cdots,\\
D_{\rm mix,implicit} &=
\hbox{coupled microscopic diffusion or selected coupled terms},\\
D_{\rm mix} &=
D_{\rm mix,explicit}+D_{\rm mix,implicit}.
\end{aligned}
\]

The exact implementation may need more than scalar \(D_{\rm mix}\) arrays
for multicomponent diffusion, but the split keeps the interface clear:
explicit/split terms keep the current path, implicit terms provide
residual and Jacobian contributions.

Any additional coupled transport source should return both its rate and
Jacobian contribution.  It should also conserve mass fraction in the
solver basis:

\[
\sum_j \mathrm{dxdt\_diff}_{j,k}=0
\quad\hbox{or}\quad
\sum_j \Phi_{j,k+1/2}=0,
\]

up to the chosen sink/projection convention.

Effort: \textbf{M} for controls and plumbing, larger if rotation's
angular-momentum coupling is pulled into the same Newton block.
\end{starbox}

# Verification Plan

\begin{verifybox}
\textbf{Finite differences are tests, not production partials.}

Raw EOS-partial check:

\[
\left.\frac{\partial f}{\partial X_j}\right|_{\rho,T}
\approx
\frac{f(X_j+\epsilon)-f(X_j-\epsilon)}{2\epsilon}.
\]

Constrained star-side check with sink species \(s\):

\[
D_j^{\rm sink}f
\approx
\frac{
 f(X_j+\epsilon,X_s-\epsilon)
 -
 f(X_j-\epsilon,X_s+\epsilon)
}{2\epsilon}.
\]

Brunt contraction check:

\[
B_{\rm comp}
=
-\frac{1}{\chi_T}
\sum_j\chi_{X_j}\frac{dX_j}{d\ln P},
\qquad
\sum_j\frac{dX_j}{d\ln P}\approx0.
\]

Run these only after explicit permission to compile/run MESA or targeted
EOS tests.
\end{verifybox}

# Implementation Order

1. Add the composition helper and source-level tests of its formulas.
2. Add Skye ideal-ion and electron `d_dxa` for `i_lnPgas` and `i_lnE`.
3. Add a Skye-local composition companion for Coulomb and phase outputs.
4. Fill all internal Skye `d_dxa(1:nv,:)` rows, but keep the public
   `num_eos_d_dxa_results = 2` until the star-side API audit is explicit.
5. Add eosDT blend `dalpha/dX_i` plumbing and alpha derivatives.
6. Revisit the Brunt branch and remove the finite-difference fallback from
   production behavior.
7. Decide which diffusion and mixing terms belong in `D_mix_explicit` and
   which belong in `D_mix_implicit`.
