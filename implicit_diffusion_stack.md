# Stack implicit diffusion on EOS composition partials

## Scope and status

Requested on 2026-10-13: back up `EbF/implicit_diffusion`, then place only
the diffusion changes on `EbF/eos_partials`. Preserve the EOS optimizations.
Do not compile or run MESA. This note is updated during the work.

- Original diffusion tip: `2b4c496d223399f52c72ddf5a5bfc38c78ed2e3f`.
- Backup: `backup/EbF-implicit-diffusion-before-eos-stack-20261013`, at the
  original diffusion tip.
- EOS parent: `2d70e611ecd1634710bbbf9fd7dfb33f60ec4e94`.
- Original diffusion upstream base: `9c454e556f91f295d8419bfd1b75ff363aacbda0`.
- EOS upstream base: `f760b022f`.
- Starting checkout: `EbF/eos_partials`, with no tracked modifications.
  Existing untracked notes, plots, model files, and outputs are user-owned.
- Source commit: `3e4838bcc`, directly on the EOS parent.
- Validation configuration commit: `40bb0a202`, directly after the source
  commit. A documentation-only commit records the completed reconstruction.
- Status: complete. `EbF/implicit_diffusion` now contains the reconstructed
  layer on `EbF/eos_partials`; this checkout is on the diffusion branch.
- The EOS branch and remote refs have not been changed. The existing
  `origin/EbF/implicit_diffusion` history differs from the rewritten local
  history. No push was performed.

## Resulting branch structure

```text
main
  EbF/eos_partials
    EbF/implicit_diffusion
```

The diffusion branch initially bundled EOS composition partials with the
implicit solver. Its first commit, `c53cf4b90`, includes both. The EOS branch
extracted the EOS layer in `e9e522ddf`, then optimized it in `2d70e611e`.
Replaying the combined commit without splitting it would duplicate that layer.

## Reconstruction invariants

1. All files under `eos/` and the EOS documentation must match the EOS parent.
   Keep its selected-row derivative evaluation and Skye basis contraction.
2. Preserve the old diffusion equations, controls, and test settings. Resolve
   shared integration files against the newer upstream code in the EOS parent.
3. Preserve normal behavior with implicit diffusion disabled. Do not introduce
   wrappers, rename routines, or move unrelated code.
4. Keep source and test configuration in the diffusion layer. Preserve the
   original notes and generated artifacts through the backup branch; do not
   duplicate them into the source change.

## Equations and code references

These are the existing diffusion equations retained by the split. The indexing
has `k=1` at the surface. `X_j` is a species mass fraction.

### Split diffusion coefficient

`star/private/implicit_Dmix.f90`, `set_Dmix_components`, retains

\[
D_{\rm mix}=D_{\rm explicit}+D_{\rm implicit},\qquad
D_{\rm explicit}\geq0,\quad D_{\rm implicit}\geq0.
\]

All three coefficients have units of cm squared per second. A full mixing pass
splits the existing total. A solver iteration holds `Dmix_explicit` fixed and
refreshes `Dmix_implicit` from the current MLT coefficient. Convection,
semiconvection, and thermohaline mixing enter the implicit bucket. The remaining
mixing contribution stays explicit.

### Composition term in the Brunt gradient

`star/private/implicit_brunt.f90`, `get_implicit_brunt_B`, evaluates at the face

\[
X_{j,f}=\alpha X_{j,k}+\beta X_{j,k-1},\qquad
\alpha=\frac{dq_{k-1}}{dq_{k-1}+dq_k},\quad \beta=1-\alpha,
\]

\[
\chi_{X_j}=\frac{P_{\rm gas}}{P_{\rm eos}}
\left(\frac{\partial\ln P_{\rm gas}}{\partial X_j}\right)_{\rho,T},\qquad
C_k=\sum_j\chi_{X_j}(X_{j,k}-X_{j,k-1}),
\]

\[
B_k=\frac{C_k-\chi_\rho\,[\ln({\rm mc}_{k-1})-\ln({\rm mc}_k)]}
{\chi_T\,[\ln P_{k-1}-\ln P_k]}.
\]

Here `mc` is the EOS mass correction. `B_1=0`. Preserve the hydrostatic pressure
difference fallback, optional Gauss path samples, and optional finite-difference
Brunt value. The existing species Jacobian holds the face `chiX` coefficients
fixed; it does not add second composition derivatives of the EOS.

The EOS-parent `get_eos_brunt_dxa_with_moments` in
`star/private/eos_support.f90` already supplies the selected rows used here.
The diffusion layer does not need to expand the stored hydro EOS derivative
array beyond its two rows.

### Mixing flux and Jacobian

`star/private/mix_info.f90`, `get_convection_sigmas`, retains

\[
\sigma_k=\frac{{\tt mix\_factor}\,D_{{\rm mix},k}
                    (4\pi r_k^2\rho_{f,k})^2}{\overline{dm}_k},\qquad
\overline{dm}_k=\frac{M_{\rm star}(dq_k+dq_{k-1})}{2}.
\]

`sigma` has units of grams per second and `sigma(1)=0`. The value uses the total
coefficient. Its structure derivatives use the implicit bucket and AD geometry;
its species derivatives use `dD_dB*d_brunt_B_dxa` multiplied by the geometric
factor. Hard clipping zeroes the coefficient derivatives. A proportional
limiter scales them with the coefficient.

`star/private/hydro_chem_eqns.f90`, `do1_chem_eqns`, preserves the coefficient
contribution to the species residual derivative:

\[
\delta R_{j,k}^{(\sigma)}=
\frac{\delta\sigma_k(X_{j,k-1}-X_{j,k})
-\delta\sigma_{k+1}(X_{j,k}-X_{j,k+1})}
{dm_k\,{\tt eqn\_scale}}.
\]

Boundary guards prevent exterior-cell references. Both structure derivatives
and off-diagonal species blocks remain present. `solver_support.f90` refreshes
the coefficients before assembling the residual.

### EOS optimization retained in the parent

`eos/private/skye.f90` retains the August 14 contraction of normalized active
number-fraction derivatives. With `Y=sum_i(X_i/A_i)` and `y_i=(X_i/A_i)/Y`, the
number-fraction contribution for an active species is

\[
\left.\frac{\partial F}{\partial X_j}\right|_{\text{number-fraction contribution}}
=\frac{F_{y_j}-\sum_i y_i F_{y_i}}{A_jY}.
\]

The composition-moment contributions are added separately. This preserves the
linear contraction and selected-row evaluation already in `2d70e611e`; no EOS
implementation files are changed by the diffusion layer.

## Integration decisions

- `turb_info.f90` passes the AD composition gradient followed by all eight
  EOS-parent thermodynamic face inputs to `do1_mlt_eval`.
- `turb_support.f90` accepts that combined interface and preserves reconstructed
  faces, pressure time centering, and sound-speed behavior. External real-valued
  MLT and semiconvection interfaces receive `%val` where the old diffusion
  implementation did so. The public `star_mlt_results` interface remains real
  valued and creates its local AD composition gradient before `Get_results`.
- `turb/private/thermohaline.f90` and `turb/public/turb.f90` retain the existing
  AD Kippenhahn, Traxler, and Brown paths. The Brown calculation differentiates
  its two implicit mode equations; it remains part of the diffusion layer.
- The newer upstream namelist readers, single-argument store routines,
  reconstructed-face profile columns/cache invalidation, and restart-aware
  rotation initialization remain present.
- Central allocation covers the new work arrays. Model/photo loading clears
  derived state and rebuilds it. Remeshing copies the split coefficients and AD
  state; full retry passes rebuild derived coefficient partials. No photo
  format change is introduced.
- The nine `implicit_diffusion_*` job controls retain defaults, declarations,
  namelist entries, and input/output assignments. Implicit diffusion remains
  off by default. Its EOS force-on and incompatible-hook checks are retained.

## Source checks

- Diffusion-only modules and solver/thermohaline integration match the backup
  exactly where the upstream code did not change.
- Reconstructed-face profile exposure, cache invalidation, and rotation restart
  changes from the EOS parent were checked in shared files.
- Authored test/work configurations match the backup apart from trailing
  whitespace and blank-line cleanup.
- All nine job controls and all nine new diffusion work fields were traced
  through their declarations, assignments, allocation, and load initialization.
- MLT/thermohaline callers were checked against the combined AD/face interface.
- EOS files, EOS documentation, and the retained EOS integration files are
  unchanged from `EbF/eos_partials`.
- Source and configuration staged whitespace checks passed. The full
  `EbF/eos_partials..EbF/implicit_diffusion` diff passes `git diff --check`.
- The EOS tip is an ancestor of the new diffusion tip; the source commit's
  direct parent is the EOS tip. The backup still points to the original
  diffusion tip.
- The reconstructed layer has no EOS implementation/doc changes or generated
  outputs. Its source and authored validation configurations are separate
  commits. The final documentation commit changes only this note.
- No compilation or MESA model execution has been performed. Numerical
  validation remains for the user, including control-off, thermohaline, Ledoux,
  `.mod` loading, photo restart, retries, and remeshing cases.

## Existing behavior to revisit separately

The original diffusion branch calls ordinary `do_brunt_B` after `set_vars` in
`finish_load_model`. This can temporarily replace the displayed implicit Brunt
value during load. It predates this split and has not been changed here.

## Independent review on 2026-10-13

Reviewed stack `0c871ca687155b431da5f6be03ddd589b839f141` against the backup
`2b4c496d223399f52c72ddf5a5bfc38c78ed2e3f` and EOS parent
`2d70e611ecd1634710bbbf9fd7dfb33f60ec4e94`. No source files were edited, and no
MESA compilation or models were run. This review updates only this note.

No lost diffusion changes or stacking regressions were found in the source
comparison. All 134 paths changed by the original branch are accounted for:

| Classification | Paths | Review result |
| --- | ---: | --- |
| Source/configuration files identical to the backup | 20 | Includes the diffusion modules, residual/Jacobian, limiting, solver refresh, remesh copies, thermohaline AD code, and six test configurations. |
| Work configurations with whitespace cleanup | 4 | Settings match the backup after ignoring trailing whitespace and final blank lines. |
| Shared source/default files integrated with the EOS parent | 13 | Original diffusion changes retained with newer face interfaces, namelist readers, profile exposure, rotation restart handling, and EOS storage. |
| Star-side files owned by the EOS parent | 6 | Retained unchanged from the EOS parent; older full-row storage and EOS-only plumbing were not reintroduced. |
| EOS implementation and documentation paths | 25 | Retained unchanged from the EOS parent. |
| EOS research/archive notes | 64 | Intentionally excluded from the diffusion layer and preserved in the backup. |
| EOS research `.gitignore` addition | 1 | Intentionally excluded; no diffusion behavior. |
| Generated `atm/test/test_output` drift | 1 | Intentionally excluded; no diffusion behavior. |

The complete `eos/` tree has the same Git object on the EOS parent and the
stack: `c4b773c3de06896c34cec0f706baa75ba57df588`. This includes file contents
and modes, and confirms that the August 14 optimizations are retained exactly.
`docs/source/eos` and the parent-owned star EOS integration files also match.

The shared MLT changes retain both the EOS parent's eight reconstructed face
inputs and the backup's AD composition-gradient input. The unique
`do1_mlt_eval` caller matches its declaration. All `Get_results` callers retain
the correct AD input; the public real-valued interface still converts locally.

All nine job controls retain complete plumbing and guards. All nine diffusion
work fields retain declaration, allocation, resizing/deallocation, and load
initialization. Load/retry passes rebuild derived state, and remesh copies are
identical to the backup. Implicit Brunt evaluates selected EOS rows into its
own full-sized scratch, so it does not require full-sized persistent hydro EOS
storage. The EOS ancestry and `git diff --check` checks passed.

### Inherited EOS diagnostic issue

`star/private/star_solver.f90:1338`, `get_grad_ad_partials`, reads
`s% d_eos_dxa(i_grad_ad,...)` when testing a composition partial. The EOS
definitions in `eos/public/eos_def.f90:108` and `:168` give `i_grad_ad=7` and
`num_eos_d_dxa_results=2`; `star/private/alloc.f90` allocates the persistent
array with those two rows. Therefore a `solver_test_partials` diagnostic for
`grad_ad` with a composition variable accesses beyond the stored rows.

This code and allocation are unchanged from `EbF/eos_partials`; the stacking
work did not introduce it. The ordinary diffusion Brunt path uses separate
scratch and is unaffected by this particular diagnostic access. A future
EOS-layer fix should request the diagnostic row into scratch rather than undo
the storage optimization. No fix is included in this review.

Source preservation is established by this review. Numerical convergence and
runtime behavior remain to be checked by the user.

## Work log

- Inspected Git history and confirmed that the two branches are siblings.
- Compared the first extracted EOS commit with the original diffusion branch;
  the EOS code was nearly identical before the August 14 optimizations.
- Created the backup branch and retained the original diffusion history.
- Created a detached temporary checkout at the EOS parent. Extracted the
  original diffusion delta relative to its upstream base and merged source
  changes onto that checkout, preserving newer upstream edits.
- Kept the EOS parent versions of `star/private/micro.f90`,
  `star/private/eos_support.f90`, `star/private/adjust_xyz.f90`,
  `star/private/init.f90`, `star/private/hydro_eqns.f90`,
  `star/private/star_solver.f90`, `star/job/run_star_support.f90`, and
  `star_data/private/star_controls_dev.inc`. Their older branch changes belong
  to the EOS layer.
- In `star/private/alloc.f90` retained two-row EOS storage and added only the
  diffusion arrays. In `star_data/public/star_data_step_work.inc` retained the
  EOS row description and added only diffusion state.
- Added the implicit-diffusion force-on behavior to the existing
  `include_eos_composition_partials` documentation in
  `star/defaults/controls_dev.defaults`.
- Preserved authored `1M_thermohaline`, `20M_pre_ms_to_core_collapse`, and
  `star/work` configurations separately from the source integration. Excluded
  EOS research notes/artifacts, archived papers, EOS-only `.gitignore`
  housekeeping, and generated `atm/test/test_output` drift from the new layer.
- Resolved the MLT interface conflicts while retaining all EOS-parent face
  inputs. Independent source review found no new integration errors in the
  solver, mixing coefficient, public API, or control paths.
- Cleaned imported configuration whitespace reported by `git diff --check`.
- Committed the diffusion-only source and initial tracking note as `3e4838bcc`,
  with `2d70e611e` as its direct parent.
- Committed the authored validation configurations as `40bb0a202` after the
  source commit. Their settings are preserved, including flags that were off
  in the original configurations.
- Updated `EbF/implicit_diffusion` with a checked old-tip comparison, then
  switched this checkout to it. Protected untracked and ignored files during
  the switch; the only existing newly tracked path was this note, whose current
  content was preserved and restored.
- Confirmed the backup tip, unchanged EOS tip, EOS ancestry, unchanged EOS
  files, and full-diff whitespace check. No MESA compilation or models were run.
- Removed the clean temporary reconstruction checkout after its commits were
  attached to the diffusion branch.

## Reviewing and continuing the stack

Review only the diffusion layer with

```sh
git diff EbF/eos_partials..EbF/implicit_diffusion
git log --oneline EbF/eos_partials..EbF/implicit_diffusion
git merge-base --is-ancestor EbF/eos_partials EbF/implicit_diffusion
git diff --exit-code EbF/eos_partials EbF/implicit_diffusion -- eos docs/source/eos
```

The first diffusion review should use `EbF/eos_partials` as its base. Further EOS
changes should be integrated beneath the diffusion layer while preserving its
diffusion-only diff. Once the EOS layer reaches main, the diffusion review can
use main as its base. The original combined implementation remains recoverable
from `backup/EbF-implicit-diffusion-before-eos-stack-20261013`.

## EOS diagnostic correction in progress

Requested on 2026-10-13: fix the `grad_ad` composition-partial diagnostic on
`EbF/eos_partials`, then update the diffusion stack. Preserve the existing
review record and the EOS optimizations. Do not compile or run MESA.

The diagnostic should request only `i_grad_ad` from `eosDT_get_dxa_rows` into
local full-sized scratch. It must not read that row from the two-row persistent
hydro storage. At fixed density and temperature, the existing species/sink
perturbation is

\[
X_p\rightarrow X_p+\epsilon,\qquad X_q\rightarrow X_q-\epsilon,
\qquad
\frac{d\nabla_{\rm ad}}{d\epsilon}
=\frac{\partial\nabla_{\rm ad}}{\partial X_p}
-\frac{\partial\nabla_{\rm ad}}{\partial X_q}.
\]

Implement this in `star/private/star_solver.f90`, `get_grad_ad_partials`, with
explicit `ierr` propagation to `test3_partials`. Keep the existing density and
temperature derivative branches and persistent EOS allocation. The diagnostic
uses the EOS parent's established blend-derivative convention documented in
`docs/source/eos/overview.rst`.

Status: fix implemented in an isolated `EbF/eos_partials` checkout. Independent
source/API review and `git diff --check` passed. Preserving this review note
before committing the EOS fix and updating the stack.

- Added a single-row EOS request into local
  `d_dxa(num_eos_basic_results,s% species)` scratch in `get_grad_ad_partials`.
- Passed current `xa`, `rho`, `T`, `lnd/ln10`, and `lnT/ln10`; the EOS API
  recomputes the composition moments from the current abundance vector.
- Retained the species-minus-sink derivative, zero neighbor contributions, and
  unchanged density/temperature diagnostic branches.
- Added explicit `ierr` propagation and a return before reading failed EOS
  output or starting numerical comparisons.
- The diagnostic requests the row regardless of
  `include_eos_composition_partials`. That control still governs production
  hydro partials. No EOS components, persistent allocation, or production
  solver equations are modified.
- The EOS API retains its existing fixed blend-weight convention and existing
  component derivative coverage; this fix does not change those conventions.
- Backed up the old EOS tip as
  `backup/EbF-eos-partials-before-grad-ad-fix-20261013` and the current stack as
  `backup/EbF-implicit-diffusion-before-grad-ad-fix-20261013`.
