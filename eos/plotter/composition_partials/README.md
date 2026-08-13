# EOS Composition Partials Validation

This is a local EOS validation plotter for the dev full-composition-partial
path.

Run it from this directory after the top-level libraries are current:

```sh
make plot
```

The Fortran driver writes the 1D finite-difference validation CSV
`data/eos_composition_partials.csv` and the 2D contour-grid CSV
`data/eos_composition_partials_contours.csv`.  The contour grid uses the
default EOS plotter domain

```text
logT    = 2..10
logRho  = -15..10
grid    = 150 x 150
```

The Python script writes contour-only PDF/PNG plots under `figures/` and a
text summary under `data/`.  The contour pages show the raw normalized
composition partials for `he4`, `c12`, and `o16`, plus a `chiX` page for the
Brunt pressure coefficient
\((P_{\rm gas}/P)(\partial\ln P_{\rm gas}/\partial X_i)_{\rho,T}\).  A
separate DFRIDR contour grid checks the constrained composition partials of
`Pgas`, `mu`, and `lnE` using `h1` as the non-plotted sink.  The no-EOS
coverage cells are masked white and are not part of the colorbar scaling.  It
also writes one combined multi-page report:

```text
figures/all_composition_partial_contours.pdf
```

The plotter removes stale `*.pdf`, `*.png`, and `*.svg` files from `figures/`
before writing the contour outputs.

## Derivative Being Checked

The EOS returns a raw fixed-\(\rho,T\) composition coefficient vector,

\[
\left.\frac{\partial q}{\partial X_i}\right|_{\rho,T}.
\]

The finite-difference check must preserve \(\sum_i X_i = 1\), so this plotter
uses `he4` as the sink species and compares

\[
\left.\frac{d q}{d X_i}\right|_{\rho,T,\sum X=1}
=
\left.\frac{\partial q}{\partial X_i}\right|_{\rho,T}
-
\left.\frac{\partial q}{\partial X_{\rm he4}}\right|_{\rho,T}
\]

against a centered Ridders extrapolation in the constrained direction

\[
X_i \rightarrow X_i + \delta X,\qquad
X_{\rm he4} \rightarrow X_{\rm he4} - \delta X .
\]

The analytic column comes from `eosDT_get_full_dxa`.  The finite-difference
value column uses the normal `eosDT_get` value path so this is a check against
the current default EOS selection, not against a special validation-only value
path.

The contour CSV separately records `raw_partial` and
`sink_projected_partial`.  The figures plot `raw_partial`.  The sink-projected
column is kept for the finite-difference validation because that is the
composition basis used by the star solver partial checks.

The DFRIDR contour CSV checks

\[
D_j^{\rm h1}q
=
\left.\frac{\partial q}{\partial X_j}\right|_{\rho,T}
-
\left.\frac{\partial q}{\partial X_{\rm h1}}\right|_{\rho,T}
\]

for \(j=\mathrm{he4},\mathrm{c12},\mathrm{o16}\).  The DFRIDR pages plot
\(\log_{10}|A-F|/\max(|A|,|F|)\), where \(A\) is the analytic constrained
partial and \(F\) is the Ridders finite-difference estimate.

For the default composition \(X=0.70\), \(Z=0.02\), the `c12` and `o16`
DFRIDR directions move the table metal coordinate through the \(Z=0.02\)
table slice:

\[
D_{\rm c12}^{\rm h1}q = D_{\rm o16}^{\rm h1}q
  = q_Z-q_X .
\]

At a table slice, MESA's analytic table derivative is the derivative of the
selected interpolation stencil, while a centered finite difference samples both
sides of the slice.  The resulting OPAL/SCVH and FreeEOS metal-direction
DFRIDR error is therefore partly a validation-stencil artifact.  The
`he4-h1` direction keeps \(Z\) fixed,

\[
D_{\rm he4}^{\rm h1}q = -q_X ,
\]

and is the cleaner check of the table \(X\) derivative at the default
composition.

## X/Z Table EOS Convention

OPAL/SCVH and FreeEOS are native table functions of hydrogen and metal
coordinates, \(q=q(X,Z)\), with helium as the remaining composition coordinate.
For the full per-isotope path, the table derivatives are mapped through the
same raw coordinates used by `chem_lib:basic_composition_info`:

\[
\frac{\partial X}{\partial X_j} = I_{\rm H}(j),
\qquad
\frac{\partial Z}{\partial X_j} = -I_{\rm H}(j)-I_{\rm He}(j),
\]

Thus

\[
\left.\frac{\partial q}{\partial X_j}\right|_{\rho,T}
=
\frac{\partial q}{\partial X}
  \frac{\partial X}{\partial X_j}
+
\frac{\partial q}{\partial Z}
  \frac{\partial Z}{\partial X_j}.
\]

For hydrogen, helium, and metal species this is respectively
\[
  q_X-q_Z,\qquad -q_Z,\qquad 0.
\]
The zero raw metal coefficient is a gauge consequence of MESA's table
coordinate definition: the table \(Z\) is the residual \(1-X_{\rm H}-X_{\rm
He}\), not the explicit sum of metal species.  Constrained derivatives
\(D_j^{\rm sink}q\) are unchanged by the equivalent row-wise constant shifts
between gauges, which is why the DFRIDR metal-direction checks above still
test \(q_Z-q_X\).

## EOS Controls

The plotter allocates a normal EOS handle with the default EOS controls.  That
is intentional for the first pass: it validates the behavior MESA will use
unless a test case supplies a custom EOS inlist.  The CSV includes the EOS
fraction rows (`frac_OPAL_SCVH`, `frac_HELM`, `frac_Skye`, `frac_PC`,
`frac_FreeEOS`, `frac_CMS`, and `frac_ideal`) so derivative failures can be
read against the active EOS component and blend regions.

## Expected Gaps

Pure PC internal composition derivatives are still deferred.  The eosDT blend
derivatives around PC are present where they have analytic alpha derivatives,
but a PC-dominated point can still show bad full-row composition derivatives
until the PC-local `MELANGE9`/mixing/phase algebra has its own analytic
composition derivative path.
