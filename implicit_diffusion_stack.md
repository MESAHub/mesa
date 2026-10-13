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
- Status: backup created; diffusion source reconstructed in
  `/private/tmp/mesa-implicit-diffusion-stack-20261013`; reviewing the shared
  MLT interfaces before updating the original branch.

## Intended branch structure

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
- Pending final staged whitespace, ancestry, and backup checks.
- No compilation or MESA model execution has been performed. Numerical
  validation remains for the user, including control-off, thermohaline, Ledoux,
  `.mod` loading, photo restart, retries, and remeshing cases.

## Existing behavior to revisit separately

The original diffusion branch calls ordinary `do_brunt_B` after `set_vars` in
`finish_load_model`. This can temporarily replace the displayed implicit Brunt
value during load. It predates this split and has not been changed here.

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
