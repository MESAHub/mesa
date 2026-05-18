# MESA SCVH low-pressure extension

This note records what the local MESA files show about the low-pressure
SCVH extension, the historical `rho = 1d-10` guard, and the local changes
needed to extend SCVH access toward `rho = 1d-15`.

## Current branch decision

Do not extend `h_tab.asc.data` or `he_tab.asc.data`.

The generated-table approach was abandoned after regenerated `eosDT`
plots still showed noisy structure in the middle OPAL/SCVH/HELM transition
region.  The current branch restores the packaged SCVH input files and
keeps the existing processed-table pressure grid:

```fortran
integer, parameter :: NlogPs = 99, NlogTs = 63
double precision, parameter :: dlogP = 0.20d0
double precision, parameter :: logP_min = -0.60d0
double precision, parameter :: logP_tail_min = -6.00d0
double precision, parameter :: logP_max = 19d0
```

The low-density extension is now analytic inside
`eos/eosDT_builder/src/scvh_core.f90`:

1. For `logP >= -0.60`, use the restored processed SCVH tables exactly.
2. For `-0.80 < logP < -0.60`, use a C1 pressure continuation that matches
   the restored bicubic table value and `d/dlogP` at `logP = -0.60`.
3. For `logP <= -0.80`, use the edge-matched ideal-gas continuation.

This keeps the table interpolation problem unchanged while still allowing
the pressure root solve to reach `rho = 1d-15`.

Older intermediate diagnostic images were removed from `notes/SCVH_eos`
after the FreeEOS fallback-row fix was verified.  The retained current
plots are listed in the latest verification section below.

The right-side `mu` structure was not caused by bicubic interpolation
inventing a feature away from the data.  Comparing the committed v51
`eosDT` table against the regenerated v53 support values shows the
feature is already present in the regenerated support table near the
high-`logQ` edge.

Two builder differences mattered:

1. The generated `logQ` grid had been moved from the committed
   `logQ_min = -10.09` grid to `logQ_min = -8.00`.  With
   `del_logQ = 0.03`, this changes the phase of every `logQ` support
   point and changes the high edge from `logQ = 5.69` to `logQ = 5.68`.
2. The old `check_results` screen rejected normal OPAL/SCVH support
   points with `gamma1 <= 1`, `gamma1 >= 2`, or `grad_ad >= 50`.  The
   regenerated table had relaxed that to reject only `gamma1 <= 0`,
   which kept finite but rough pressure-ionization values on the right
   side.  In the approximate right-side box
   `0 <= log10(rho) <= 3`, `3 <= log10(T) <= 4.6`, 1202 regenerated
   support cells would have been rejected by the old screen.  None of
   those are in the low-pressure analytic-tail exemption.

The current code keeps the old derivative screen for normal OPAL/SCVH
support points and exempts only the analytic low-pressure tail:

```fortran
lowP_extension = logPgas_opalscvh < -0.6d0

bad = gamma1_opalscvh <= 0d0 .or. &
      (.not. lowP_extension .and. &
       (gamma1_opalscvh <= 1d0 .or. gamma1_opalscvh >= 2d0 .or. grad_ad_opalscvh >= 50d0))
```

For the low-density diagnostic box
`-15 <= log10(rho) <= -8`, `2.1 <= log10(T) <= 3.3`, 840 of the
old-screen rejections are below `logPgas = -0.6` and remain on the
analytic-tail branch; 177 still fail outside that tail and should fall
back as before.

Older guarded-v53 plotter images were removed from `notes/SCVH_eos`.

The plotter rejected the stale binary caches with `num_logQs` and
`eos_logQ_min` mismatches and rewrote the v53 caches from the rebuilt
text tables.  The rebuilt table header for `mesa-eosDT_02z60x.data` is:

```text
version = 53
num logQs = 527
logQ min = -10.0900
logQ max = 5.6900
```

The rebuilt v53 support grid is now on the committed v51 `logQ` grid.
It is not bitwise-equivalent to the v51 data because the strict screen
now routes many pressure-ionization support points through the current
HELM/fallback path.  For `mesa-eosDT_02z60x.data`, compared cell-by-cell
with the committed v51 table:

```text
right stripe, 0 <= log10(rho) <= 3 and 3 <= log10(T) <= 4.6:
  cells = 3577
  |delta mu| max/95%/median = 1.28894 / 0.88735 / 0.04623
  rebuilt HELM nonzero/max = 1210 / 1.0

right box, -0.5 <= log10(rho) <= 4.5 and 3 <= log10(T) <= 4.8:
  cells = 6097
  |delta mu| max/95%/median = 1.77469 / 1.53301 / 0.03629
  rebuilt HELM nonzero/max = 2264 / 1.0
```

So the grid-phase problem is fixed, but the rebuilt table still differs
from the committed table in screened pressure-ionization cells because
the current fallback EOS path is not the same as the historical v51
generator path.

Follow-up: the large right-side differences are specifically cells where
the guarded rebuild used current HELM with `mu = 0.66835`, while the
committed v51 support table has molecular/neutral OPAL/SCVH-like values
around `mu = 1.9..2.3`.  The old table therefore either did not fall
back there, or did not fall back to today's complete-ionization HELM path.
The builder policy is now:

1. Keep finite OPAL/SCVH pressure-ionization support points, even when
   `gamma1` or `grad_ad` is large.
2. Use the SCVH replacement for OPAL failures only on the low-density
   side, `logRho < logRho5`, so high-density/right-edge OPAL failures do
   not incorrectly keep SCVH.
3. Restore the old SCVH edge retry condition to `info == -2`, rather
   than retrying all negative SCVH statuses.

The generated `eosDT` file version was set to 52 because this change is
intended as a single replacement for the packaged version 51 tables.

Older no-finite-fallback plotter images were removed from `notes/SCVH_eos`.

The plotter rejected the previous binary caches with
`read cache failed for version_in`.  For `mesa-eosDT_02z60x.data`, the
version 52 table comparison against committed v51
is:

```text
right stripe, 0 <= log10(rho) <= 3 and 3 <= log10(T) <= 4.6:
  cells = 3577
  |delta mu| max/95%/median = 1.22261 / 0.170812 / 0.0
  rebuilt HELM nonzero/max = 5 / 1.0

right narrow visual, -0.2 <= log10(rho) <= 1.5 and 3.4 <= log10(T) <= 4.8:
  cells = 3610
  |delta mu| max/95%/median = 0.38484 / 0.135336 / 0.0
  rebuilt HELM nonzero/max = 0 / 0.0
```

This confirms the right-side artifact was mainly due to finite
pressure-ionization cells being routed to current HELM.

The plot uses one-sided table-edge slopes for visualization.  The Fortran
implementation uses the bicubic edge derivatives already returned by
`interp_evbicub_db`.

Implementation note: the analytic branch reuses the existing
`interp_value` path at `logP_min`, then applies the tail formula.  There
is no second bicubic interpolation path and no generated-table helper
script.

Builder region note: `eos_regions_code.dek` uses `logT < logT7`, not
`logT <= logT7`, so the first table row at `logT = 2.10` belongs to the
OPAL/SCVH side.  With the old inclusive test, the regenerated builder
sent `logT = 2.10`, `logrho = -15` to HELM, which is below HELM's valid
range and stopped in `helmeos2`.

### Analytic tail equations

Let

```text
p  = log10(P)
p0 = -0.60
pd = -0.80
L  = p0 - pd
x  = (p - pd)/L
```

At the processed SCVH edge, evaluate the restored bicubic table:

```text
y0 = y(p0, logT)
m0 = (dy/dp)(p0, logT)
```

For density, internal energy, and composition, the edge-matched ideal
continuations are:

```text
logrho_id(p) = logrho0 + p - p0
logu_id(p)   = logu0
comp_id(p)   = comp0
```

For the ideal-tail density derivative, use
`dlogrho_id/dlogT = -1`.  The edge-derived
`Rgas = P0/(rho0 T)` is a species constant in the ideal tail, so it is
not differentiated through the bicubic edge derivative.

For entropy:

```text
Rgas = 10^(p0 - logrho0 - logT)
S_id(p) = S0 + Rgas ln(10) (p0 - p)
logs_id(p) = log10(S_id(p))
```

The C1 pressure correction is:

```text
g(x) = x^3 - x^2
M    = (m0 - mid0) L
y(p) = y_id(p) + M g(x),        pd < p < p0
y(p) = y_id(p),                 p <= pd
```

Here `mid0` is the ideal derivative at `p0`: `1` for `logrho`, `0` for
`logu` and frozen composition, and `-Rgas/S0` for `logs`.  Since
`g(0) = g(1) = 0`, `g'(0) = 0`, and `g'(1) = 1`, the tail is value- and
pressure-derivative-continuous at both `pd` and `p0`.  This is C1 in
pressure, not a C2 reconstruction.

## Short version

The original SCVH H/He tables are pressure-temperature tables, not
density-temperature tables:

```text
logT, logP -> rho, S, U, composition fractions
```

In the first SCVH-derived files imported into MESA, the cold
low-pressure edge is at about `logP = 4`; the `logP = -0.60` tail was
not present at first import.  Before this local extension, later MESA
releases shipped processed SCVH files in which the cold low-pressure
tail had already been extended down to `logP = -0.60`.  The coldest
values in that added region are consistent with an ideal neutral gas
continuation:

Thus the local archive evidence says the table was extended after SCVH
was introduced into MESA, not before; the extension appeared somewhere
between the r119 import and the r2025 input-data text-base state.

```text
H:  molecular H2, diatomic ideal gas
He: neutral He, monatomic ideal gas
```

Before this local extension, the `rho = 1d-10` guard in `scvh_core.f90`
was a MESA runtime safety guard for the SCVH density-to-pressure
inversion.  It is not an original SCVH density boundary.  The MESA
repository does not show a note proving whether the guard motivated the
`logP = -0.60` processed-table edge, or whether the two choices were made
independently.  What the files show is only that the choices are
numerically consistent.  In this branch the processed SCVH input files are
kept unchanged, while `scvh_core.f90` adds an analytic tail and lowers the
runtime guard to `Rho_min = 1d-15`.

## Relevant local files

Raw/public-SCVH-like inputs:

```text
eos/eosDT_builder/eos_input_data/scvh/h_tab_initial.dat
eos/eosDT_builder/eos_input_data/scvh/he_tab_initial.dat
```

MESA-processed SCVH inputs:

```text
eos/eosDT_builder/eos_input_data/scvh/h_tab.dat
eos/eosDT_builder/eos_input_data/scvh/he_tab.dat
```

Final 5-column files read by current `scvh_core.f90`:

```text
eos/eosDT_builder/eos_input_data/scvh/h_tab.asc.data
eos/eosDT_builder/eos_input_data/scvh/he_tab.asc.data
```

Code that reads the processed SCVH files:

```text
eos/eosDT_builder/src/scvh_core.f90
```

The pressure grid remains hard-coded as:

```fortran
integer, parameter :: NlogPs = 99, NlogTs = 63
double precision, parameter :: dlogP = 0.20d0
double precision, parameter :: logP_min = -0.60d0
double precision, parameter :: logP_tail_min = -6.00d0
double precision, parameter :: logP_max = 19d0
double precision, parameter :: dlogT = 0.08d0
double precision, parameter :: logT_min = 2.10d0
double precision, parameter :: logT_max = 7.06d0
```

The processed SCVH input grid-size identity is:

```text
NlogPs = (logP_max - logP_min)/dlogP + 1
       = (19 - (-0.60))/0.20 + 1
       = 99
```

The analytic `logP_tail_min = -6.00` branch is not an added table grid;
it is evaluated in `scvh_core.f90` when the requested pressure falls
below the processed edge.

The runtime density guard is:

```fortran
parameter(Rho_min=1.0d-15)

if (Rho < Rho_min) then
   Rho = Rho_min
   logRho = log10(Rho)
end if
```

## What is original SCVH and what is MESA processing?

The local `h_tab_initial.dat` starts at `logP = 4.00` for the coldest
temperature row:

```text
logT = 2.10, logP = 4.00, logrho = -5.7154
```

The local `he_tab_initial.dat` likewise starts at `logP = 4.00`:

```text
logT = 2.10, logP = 4.00, logrho = -5.4175
```

So the cold low-pressure region below `logP = 4` is not present in the
initial SCVH-derived files.  In particular, the SCVH tables were not
already extended to `logP = -0.60` when SCVH first entered MESA in SVN
r119.

The processed MESA files extend below that.  In the superseded
generated-table attempt, the old 11-column intermediate files were left
unchanged, while the final `.asc.data` files read by `scvh_core.f90`
were extended farther:

```text
h_tab.dat      extends down to logP = -0.60
he_tab.dat     extends down to logP = -0.60
h_tab.asc.data extends down to logP = -6.00  (superseded)
he_tab.asc.data extends down to logP = -6.00 (superseded)
```

The current branch restores `h_tab.asc.data` and `he_tab.asc.data` to the
packaged `logP = -0.60` edge and implements the lower-pressure tail in
`scvh_core.f90` instead.

The current source reads the already-processed `.asc.data` files.  The
historical script or routine that produced the `logP = 4` to `-0.60`
extension from the initial files is not present as an obvious standalone
reproducer in the current tree.

The local checked-in `.svn` metadata under:

```text
eos/eosDT_builder/eos_input_data/scvh/.svn
```

records an intermediate SVN working-copy state for the input data:

```text
Last changed rev: 2025
Last changed date: 2010-01-06T01:27:30.336251Z
Last changed author: aaron_dotter
URL: https://mesa.svn.sourceforge.net/svnroot/mesa/trunk/eos/eosDT_builder/eos_input_data/scvh
```

The `.svn/text-base` copies match the local `h_tab_initial.dat`,
`he_tab_initial.dat`, `h_tab.dat`, and `he_tab.dat` files.  That proves
the `*_initial.dat` raw tables and the already-extended `*.dat` tables
coexisted in the MESA SVN input-data tree by r2025.  It does not prove
that the extension existed before SCVH was introduced to MESA; the
fetched r119 import shows the opposite for the first imported SCVH
files.

There appear to be two historical preprocessing stages:

1. `h_tab_initial.dat` / `he_tab_initial.dat` to `h_tab.dat` /
   `he_tab.dat`.

   This is where the cold low-pressure rows below `logP = 4` were added.
   For the low-temperature rows, `h_tab.dat` has exactly 23 more pressure
   points than `h_tab_initial.dat`, corresponding to:

   ```text
   logP = -0.60, -0.40, ..., 3.80
   ```

   The local r2025 `.svn` metadata shows this stage had already happened
   by 2010-01-06.  The current MESA tree and the local fetched
   SCVH-grep commit archive do not contain an obvious routine that reads
   `h_tab_initial.dat`, adds these 23 low-pressure rows, merges the
   pressure ionization pieces, and writes `h_tab.dat`.

2. `h_tab.dat` / `he_tab.dat` to `h_tab.asc.data` / `he_tab.asc.data`.

   This stage is still partly represented in `scvh_core.f90`.  The
   current `read_file_for_scvh` routine can read the first five value
   columns from files shaped like `h_tab.dat`, `fill_and_smooth` fills
   missing rectangular-grid regions, and the disabled
   `write_file_for_scvh` path can write `.asc.data`-style output.  The
   source currently reads `.asc.data` directly, so this writer is not part
   of normal operation.

Thus, the exact first-stage historical generator for the `logP = 4` to
`-0.60` tail is not available as a current MESA command to rerun.  The
local files prove the result of that step and bound its history between
the original r119 import and the r2025 input-data state, but they do not
show the formula-bearing routine.

## Check of older local MESA versions

The oldest local MESA release tree under `/Users/owner/Documents/Software`
is:

```text
/Users/owner/Documents/Software/mesa-r11701
```

That tree does not contain the missing first-stage preprocessing routine.
Before the local `logP_min = -6` extension, its SCVH input tarball was
byte-identical to the then-current one:

```text
MD5 (mesa-r11701/eos/eosDT_builder/eos_input_data.tar.xz) = 029396c7fc3ece7496a819dfd4e8deda
MD5 (pre-extension/eos/eosDT_builder/eos_input_data.tar.xz) = 029396c7fc3ece7496a819dfd4e8deda
```

So `mesa-r11701` already shipped the same processed SCVH input data,
including the same `logP = -0.60` low-pressure extension.

The old PT-side SCVH reader in:

```text
/Users/owner/Documents/Software/mesa-r11701/eos/eosPT_builder/src/scvh_eval.f90
```

also assumes the processed table.  It hard-codes:

```fortran
double precision, parameter :: max_scvh_logP = 19d0, min_scvh_logP = -0.6d0
integer, parameter :: nt=63, np=99
```

and reads:

```fortran
scvh/h_tab.dat
scvh/he_tab.dat
```

with the pressure grid:

```fortran
plog(j) = -0.60d0 + (j-1)*0.20d0
```

This is a reader/interpolator for the already-extended table, not a
generator for `h_tab.dat` from `h_tab_initial.dat`.

Searches across the local older release trees:

```text
mesa-r11701
mesa-r12115
mesa-r12778
mesa-r15140
mesa-main
mesa-r26.04.1
```

found no source references to:

```text
h_tab_initial
he_tab_initial
h_tab_p1
h_tab_p2
```

Thus the low-pressure-tail generator is not present in the older local
MESA releases checked so far.  It likely predates the packaged release
history present locally, lived outside the released MESA source tree, or
was an old preprocessing artifact not retained with the builder code.

## SVN mirror findings

SCVH-related commits were fetched from the old cgit/SVN mirror into:

```text
notes/old_mesa_eos_info/svn_scvh_commits
```

The curated local index is:

```text
notes/old_mesa_eos_info/svn_scvh_commits/important_scvh_commits.md
```

Two fetched commits, plus one local `.svn` working-copy data state, are
especially relevant to the low-pressure-extension question:

1. `0a3649bd8c55b5f164278139a68372c04c6850ae`

   ```text
   2007-01-29, SVN r119, add scvh data
   ```

   This imports:

   ```text
   eos/eos_builder/eos_input_data/scvh/h_tab.dat
   eos/eos_builder/eos_input_data/scvh/h_tab_p1.dat
   eos/eos_builder/eos_input_data/scvh/h_tab_p2.dat
   eos/eos_builder/eos_input_data/scvh/he_tab.dat
   ```

   The imported H table starts:

   ```text
   logT = 2.10, nps = 30
   first pressure point: logP = 4.00
   ```

   So the original imported `h_tab.dat` in the fetched history did not
   yet contain the `logP = -0.60` tail.

2. Local `.svn` working-copy metadata for:

   ```text
   eos/eosDT_builder/eos_input_data/scvh
   ```

   records:

   ```text
   SVN r2025, 2010-01-06, author aaron_dotter
   ```

   In that local SVN text-base state, the raw/public-SCVH-like files:

   ```text
   h_tab_initial.dat
   he_tab_initial.dat
   ```

   start at `logP = 4.00`, while the processed files:

   ```text
   h_tab.dat
   he_tab.dat
   ```

   already start at `logP = -0.60` for the cold rows.  The local
   `svn_scvh_commits` archive does not contain the r2025 patch, so this
   identifies the data state but not the generator routine.

3. `ece7a79430f66e84687d5b0ce52288a211b7dc84`

   ```text
   2010-04-15, SVN r2260, turn on entropy of mixing in scvh eos
   ```

   This creates:

   ```text
   eos/eosDT_builder/h_tab.asc.data
   eos/eosDT_builder/he_tab.asc.data
   ```

   The new H `.asc.data` file starts:

   ```text
   logT = 2.10, nps = 99
   first pressure point: logP = -0.60
   ```

   That is the first fetched commit where the processed low-pressure
   SCVH `.asc.data` tail appears explicitly.  The same patch changes
   `make_SCVH_plots.f` to use `logP_min = -0.6d0` in the SCVH plotting
   ranges and changes `scvh_core.f` to read `.asc.data` directly.

The local history therefore says:

```text
r119:  SCVH enters MESA without the -0.60 low-pressure tail.
r2025: local SVN text-base has both initial and already-extended .dat files.
r2260: fetched patch adds generated .asc.data files with the -0.60 tail.
```

The r2025 local data state narrows the historical source of the processed
tail, but it still does not provide a clean standalone recipe for
generating the `logP = 4.00` to `-0.60` extension from the 2007 data
files.  The useful missing historical target is the r2025 patch or the
full source state immediately around r2025.

### What this implies for the superseded local `-6.00` table implementation

The local tables do let us test the simple low-pressure continuation
against the historical `*_initial.dat` to `*.dat` data change.  Comparing
the added `logP = -0.60, -0.40, ..., 3.80` rows in `h_tab.dat` and
`he_tab.dat` to the original `logP = 4.00` edge shows:

```text
For cold neutral rows:
  logrho(logP) = logrho(4.00) + logP - 4.00
  logu(logP)   = logu(4.00)
  S(logP)      = S(4.00) + R_species ln(P(4.00)/P)
  composition  = composition(4.00)
```

For those rows the mismatch from this model is only roundoff/table
formatting.  At `logT = 2.10`, for example, the maximum `logrho` mismatch
is below `2d-15`, `logu` and the composition columns are unchanged, and
the entropy mismatch is about `1d-4` in relative `S`.

That simple continuation is not a faithful description of the whole
historical `logP = 4.00` to `-0.60` extension.  In the H dissociation /
ionization region, the local table comparison shows changes in
composition and internal energy across the added low-pressure points.
For example, around `logT = 3.30` in H the added rows differ from a
frozen-edge continuation by order-unity amounts in the composition
columns and by nearly one dex in `logu`.  Around `logT = 4.02` the H and
He rows also differ substantially from a frozen-edge continuation.

The superseded `extend_scvh_low_p_tail.py` implementation was therefore
historically supported only as a cold-neutral ideal-gas boundary
continuation below the already-processed `logP = -0.60` edge.  It should
not be described as recovering the original MESA preprocessing routine
for all of the `logP = 4.00` to `-0.60` region.  If we need correctness
through the dissociation / ionization rows below `logP = -0.60`, the
missing item is still the r2025 generator or an independently justified
chemical equilibrium continuation.

### Provenance of the historical preprocessing hook

The closest historical source for the preprocessing path we reconstructed
is the same 2010-04-15 commit:

```text
Commit:  ece7a79430f66e84687d5b0ce52288a211b7dc84
SVN rev: 2260
Date:    2010-04-15
Subject: turn on entropy of mixing in scvh eos
Author:  bill-paxton
```

In that patch, the original Fortran file name was:

```text
eos/eosDT_builder/src/make_SCVH_plots.f
```

and the relevant entry printed:

```text
make_SCVH_files
```

The same patch added the generated processed value files:

```text
eos/eosDT_builder/h_tab.asc.data
eos/eosDT_builder/he_tab.asc.data
```

and added a disabled writer hook in:

```text
eos/eosDT_builder/src/scvh_core.f
```

with the subroutine name:

```text
write_file_for_scvh
```

writing:

```text
scvh/h_tab.new.data
scvh/he_tab.new.data
```

That is the original-name/provenance trail for the processed `.asc.data`
path.  The local Python helper should be described as sourced from that
historical MESA path and from the ideal-gas low-pressure continuation
visible in the generated r2260 `.asc.data` files.  It is not known to be
the exact historical routine, because the standalone first-stage
generator that made the `logP = 4.00` to `-0.60` tail has not been found.
The local Python helper:

```text
eos/eosDT_builder/extend_scvh_low_p_tail.py
```

is not a recovered SVN script.  It is a reconstruction of the missing
first-stage low-pressure-tail generation, written to extend the already
processed `.asc.data` files using the ideal-gas continuation visible in
the 2010 generated data.

## Entropy-of-mixing history

The entropy-of-mixing path was checked separately from the low-pressure
tail generation because it answers a different question: whether the
active H/He SCVH mixing implementation was later removed.

In the active DT-builder path, it was not removed.  Current
`eos/eosDT_builder/src/scvh_core.f90` still calls
`entropy_of_mixing`, adds `smix` to the mixed entropy, and includes
`d_smix_dT` and `d_smix_dP` in the entropy derivatives:

```text
call entropy_of_mixing(...)
!smix = 0; d_smix_dT = 0; d_smix_dP = 0
entr = xmassh1*entr_h + xmasshe4*entr_he + smix
dsdt_cp_hhe = (... + T*d_smix_dT)/entr
dsdpress_ct_hhe = (... + P*d_smix_dP)/entr
```

The relevant code references are:

```text
eos/eosDT_builder/src/scvh_core.f90:794
eos/eosDT_builder/src/scvh_core.f90:798
eos/eosDT_builder/src/scvh_core.f90:801
eos/eosDT_builder/src/scvh_core.f90:816
eos/eosDT_builder/src/scvh_core.f90:817
eos/eosDT_builder/src/scvh_core.f90:1064
```

The fetched 2010-04-15 SVN patch (`SVN r2260`, subject `turn on entropy
of mixing in scvh eos`) is where this DT-builder contribution was turned
on.  The older DT code computed `smix` and its derivatives and then
zeroed them with the comment:

```text
SKIP SMIX FOR NOW.  IT MESSES UP GRADAD.
```

The patch comments out that zeroing line instead:

```fortran
!smix = 0; d_smix_dT = 0; d_smix_dP = 0
```

The later deletion found in the fetched history is not a deletion of the
active `entropy_of_mixing` routine.  The 2020-05-14 SVN patch
(`SVN r13653`, subject `delete SCVH_DT stuff. we will have to wait for
CMS.`) deletes the standalone SCVH_DT data product/script and comments
out `scvh_only = .true.` in the DT builder.  It does not remove the
mixed OPAL/SCVH `eosDT_builder/src/scvh_core` entropy-of-mixing code.

The older PT-side evaluator is different.  Old local releases include
`eosPT_builder/src/scvh_eval.f90`, and the 2020 SCVH_PT import includes
`eosSCVH_PT_builder/src/scvh_eval.f90`.  Those files compute an
approximate `smix`, but do not add it to `entr`; comments say the
partials are missing and that FXT said to skip `smix` because the
partials would require differencing the tables.  That skipped PT-side
implementation is not the active DT-builder path used to generate the
packaged OPAL/SCVH/eosDT tables.

Local checks across available old releases (`r11701`, `r12115`,
`r12778`, `r15140`, `mesa-main`, and `r26.04.1`) show the active
`eosDT_builder/src/scvh_core` entropy-of-mixing path still present.
Current git history only reports the initial MESA import and the later
`.f` to `.f90` conversion for searches around `entropy_of_mixing` or
`smix = 0`; it does not show a later removal from the active DT-builder
source.

## Superseded: what would be needed to extend the table farther

To extend the processed SCVH input in the same spirit as the existing
`logP = 4.00` to `-0.60` tail, the work had to be done at three levels:

1. Generate a new processed H/He pressure tail.
2. Teach `scvh_core.f90` about the larger pressure grid and lower
   density inversion guard.
3. Regenerate the packaged `eosDT` tables and move the runtime
   OPAL/SCVH selection limits so the new region is actually used.

Before this implementation, the SCVH reader assumed:

```fortran
integer, parameter :: NlogPs = 99, NlogTs = 63
double precision, parameter :: dlogP = 0.20d0
double precision, parameter :: logP_min = -0.60d0
double precision, parameter :: logP_max = 19d0
```

with:

```text
NlogPs = (logP_max - logP_min)/dlogP + 1
```

For example, extending to `logP_min = -6.00` at the same spacing would
require:

```text
NlogPs = (19 - (-6))/0.20 + 1 = 126
```

The low-pressure generated values should preserve the existing
processed-table composition convention at the old low-pressure edge for
each fixed `logT`.  In the coldest H row this is:

```text
H table columns at logT = 2.10, logP = -0.60:
  0.00000E+00  1.00000E+00

He table columns at logT = 2.10, logP = -0.60:
  1.00000E+00  0.00000E+00
```

The extension keeps those column values fixed below `logP = -0.60`
rather than changing the historical column convention.  The important
closure condition for avoiding spurious ionization is that the neutral
plus molecular/ion fractions remain inherited from the processed edge:

```text
for the cold neutral rows:
x_H+   = 1 - (x_H2 + x_H) = 0
x_He++ = 1 - (x_He+ + x_He) = 0
```

The density entries should follow the ideal neutral-gas relations used
by the existing low-pressure tail:

```text
H2: logrho = logP - logT - log10(k_B/(2 m_H))
He: logrho = logP - logT - log10(k_B/m_He)
```

Numerically, from the existing table:

```text
log10(k_B/(2 m_H)) ~= 7.6154
log10(k_B/m_He)    ~= 7.3175
```

Thus `rho = 1d-15` at the cold edge requires roughly:

```text
H2, logT = 2.10: logP = -15 + 2.10 + 7.6154 = -5.2846
He, logT = 2.10: logP = -15 + 2.10 + 7.3175 = -5.5825
```

So a `logP_min` near `-6` would cover `rho = 1d-15` at `logT = 2.10`
for both pure-H and pure-He limits.

For the energy and entropy entries, the least invasive reconstruction is
to splice onto the pre-extension processed edge at each fixed `logT`:

```text
u(logP < logP0, T) = u(logP0, T)
S(logP < logP0, T) = S(logP0, T) - R_species ln(P/P0)
```

where `logP0` was the first existing processed pressure point, `-0.60`
for low-temperature rows, and:

```text
R_H2 = k_B/(2 m_H)
R_He = k_B/m_He
```

The table stores `log10(S)` and `log10(u)`, so this splice has to be
done in linear `S` and `u`, then converted back to logarithms.  This
keeps the new tail continuous at the old edge and gives the ideal-gas
pressure derivative for entropy.  A more ambitious version could
reconstruct the full original first-stage preprocessing from
`h_tab_initial.dat` / `he_tab_initial.dat`, but that routine has not yet
been found in the current or local old MESA trees.

The source changes implied by a farther pressure tail are:

```text
eos/eosDT_builder/src/scvh_core.f90
  NlogPs
  logP_min
  Rho_min
  write_file_for_scvh hard-coded nps values, if the writer is reused
```

For `logP_min = -6`, a `Rho_min` of `1d-15` is consistent with the
table; the pre-extension `Rho_min = 1d-10` would still have clamped away
the new low-density region before the pressure solve.

The generated packaged table also had to sample the lower density:

```text
eos/eosDT_builder/src/create_eos_files.f90
  logRho_min = -14d0
```

needed to be lowered if the final `eosDT` product was meant to include
`logRho = -15`.

Finally, runtime selection limits also matter.  Before this
implementation, the OPAL/SCVH selection logic had low-density and
low-`logQ` cutoffs:

```text
eos/defaults/eos.defaults
  logRho_min_OPAL_SCVH_limit = -14.299d0
  logQ_min_OPAL_SCVH = -8.0d0

eos/private/eosdt_eval.f90
  logRho5 = rq% logRho_min_OPAL_SCVH_limit
  logRho6 = logRho5 - 1d-3
  logRho7 = -14.90d0
  logRho8 = -14.99d0
```

At `logT = 3` and `logRho = -15`, `logQ = -9`, so the pre-extension
`logQ_min_OPAL_SCVH = -8` would have rejected OPAL/SCVH there even if
the pressure table had been extended.  Lowering
`logRho_min_OPAL_SCVH_limit` also requires moving the hard-coded
`logRho7` and `logRho8` values so the required monotonic ordering of the
region boundaries is preserved.

In short, extending the SCVH tail is feasible, but it is not just one
constant.  The minimum consistent implementation is:

```text
new h_tab.asc.data / he_tab.asc.data tail
new scvh_core pressure-grid constants
lower scvh_core Rho_min
lower eosDT builder logRho_min, if the packaged table must reach there
move OPAL/SCVH runtime logRho and logQ cutoffs
regenerate eosDT_data.tar.xz and invalidate old caches
```

## Superseded implemented low-pressure table extension

The working implementation targets `rho = 1d-15` at the low-temperature
ideal-gas boundary by extending the processed SCVH pressure tail to:

```text
logP_min = -6.00
dlogP = 0.20
NlogPs = 126
```

The new local preprocessing helper is:

```text
eos/eosDT_builder/extend_scvh_low_p_tail.py
```

It reads the current processed files:

```text
eos/eosDT_builder/eos_input_data/scvh/h_tab.asc.data
eos/eosDT_builder/eos_input_data/scvh/he_tab.asc.data
```

and prepends ideal-gas entries below the previous `logP = -0.60` edge.
This helper is sourced from the historical r2260 processed-SCVH path and
the ideal-gas tail visible in those generated data files.  It is a
reconstructed replacement for the missing first-stage preprocessing step,
not a byte-for-byte recovered SVN script.  The SVN provenance for the
historical path is:

```text
SVN r2260, 2010-04-15, ece7a79430f66e84687d5b0ce52288a211b7dc84
  original nearby driver: eos/eosDT_builder/src/make_SCVH_plots.f
  printed entry:          make_SCVH_files
  writer hook:            write_file_for_scvh in scvh_core.f
  generated files:        h_tab.asc.data, he_tab.asc.data
```

### Continuation formula used by the helper

The `.asc.data` table stores base-10 logarithms for pressure,
temperature, density, entropy, and internal energy:

```text
ell_P   = log10(P)
ell_T   = log10(T)
ell_rho = log10(rho)
ell_S   = log10(S)
ell_u   = log10(u)
```

For each fixed-temperature row, let the first existing processed point be
the splice point:

```text
ell_P0   = first.logp
ell_T    = row.logt
ell_rho0 = first.logrho
ell_S0   = first.logs
ell_u0   = first.logu
c1_0     = first.comp1
c2_0     = first.comp2
```

For the current low-temperature rows, `ell_P0 = -0.60`.  The script only
extends rows whose first point is already below the original SCVH-like
`logP = 4` edge.  Rows that still start at `logP ~= 4` are left for the
existing in-memory `fill_and_smooth` path in `scvh_core.f90`.

The ideal-gas constant per gram is inferred from the existing edge point:

```text
R_species = P0/(rho0 T)
          = 10^(ell_P0 - ell_rho0 - ell_T)
```

This is the same relation as:

```text
P = rho R_species T
```

so holding `R_species` fixed gives the density continuation:

```text
ell_rho(ell_P) = ell_rho0 + ell_P - ell_P0
```

The ideal reference tail holds the energy fixed at fixed temperature:

```text
u(ell_P)     = u0
ell_u(ell_P) = ell_u0
```

The entropy is continued in linear entropy, not in `log10(S)`:

```text
S0 = 10^ell_S0

S(ell_P) = S0 + R_species ln(P0/P)
         = S0 + R_species ln(10) (ell_P0 - ell_P)

ell_S(ell_P) = log10(S(ell_P))
```

The ideal reference tail also holds the composition columns fixed:

```text
c1_ideal(ell_P) = c1_0
c2_ideal(ell_P) = c2_0
```

In the cold H rows this preserves the processed molecular/neutral
fractions.  In the cold He rows it preserves neutral helium.  The closure
relations therefore remain neutral in the artificial tail:

```text
x_H+   = 1 - (x_H2 + x_H) = 0
x_He++ = 1 - (x_He+ + x_He) = 0
```

The first version used these ideal values directly.  That made the tail
continuous in value at `ell_P0`, but it forced zero pressure derivative
in `c1`, `c2`, and `ell_u` below the splice, even when the already
processed table above `ell_P0` had a nonzero slope.  Because
`scvh_core.f90` builds bicubic interpolation data from the full table,
that value-only splice is not a harmless local detail.

The current helper therefore uses the ideal tail as a far-tail reference
and adds a smooth first-derivative correction near the splice.  For each
stored quantity

```text
y in {c1, c2, ell_rho, ell_S, ell_u}
```

let:

```text
m_old = [y(ell_P0 + dlogP) - y(ell_P0)]/dlogP
m_id  = dy_ideal/dell_P at ell_P0
```

with:

```text
m_id(c1)      = 0
m_id(c2)      = 0
m_id(ell_rho) = 1
m_id(ell_u)   = 0
m_id(ell_S)   = -R_species/S0
```

To keep the actual density floor ideal, the correction is not allowed to
persist all the way down to `ell_Pmin` when the row reaches
`log10(rho) = -15` at a higher pressure.  Define the ideal-pressure
anchor:

```text
ell_rho_anchor = -15.00
ell_P_id = ell_P0 + (ell_rho_anchor - ell_rho0)
ell_P_decay = min(ell_P0 - dlogP, max(ell_Pmin, ell_P_id))
```

For `ell_P_decay < ell_P < ell_P0`, define:

```text
q = (ell_P0 - ell_P)/(ell_P0 - ell_P_decay)
w(q) = 1 - (10 q^3 - 15 q^4 + 6 q^5)
```

The helper also enforces discrete one-cell slope continuity at the first
new grid point:

```text
ell_P1 = ell_P0 - dlogP
y_match = y(ell_P0) - m_old dlogP
```

This makes:

```text
[y(ell_P0) - y(ell_P1)]/dlogP
   = [y(ell_P0 + dlogP) - y(ell_P0)]/dlogP
```

to table precision.  The extra adjustment is an endpoint bump:

```text
t = (ell_P - ell_P_decay)/(ell_P0 - ell_P_decay)
b(t) = t^2 (1 - t)^2
```

with an amplitude chosen so that `y(ell_P1) = y_match`.  In code form:

```text
y(ell_P) = y_ideal(ell_P)
         + (m_old - m_id) (ell_P - ell_P0) w(q)
         + A b(t)
```

where `A` is solved independently for each table column and row.  Since
`b` and `db/dell_P` vanish at both endpoints, this bump does not change
the value or derivative at the splice, and it does not change the ideal
tail at `ell_P_decay`.

For `ell_P <= ell_P_decay`, the helper writes the ideal-reference value
directly.

This preserves the edge value, matches the one-sided pressure derivative
from the old processed table at `ell_P0`, matches the first discrete
finite-difference slope used by the sampled table, and decays back to the
ideal tail by the row's `rho = 1d-15` anchor, or by `ell_Pmin` if that is
lower.  The composition pair is then limited to keep
`0 <= c1`, `0 <= c2`, and `c1 + c2 <= 1`.

The pressure points prepended by the helper are:

```text
ell_P = logP_min, logP_min + dlogP, ..., ell_P0 - dlogP
```

with:

```text
logP_min = -6.00
dlogP = 0.20
```

The existing `ell_P0` point is retained from the input file, so the new
tail is continuous at the splice point.

Inside the far tail, where the correction has decayed away, the analytic
pressure derivatives return to the ideal-reference values:

```text
d ell_rho / d ell_P = 1
d ell_u   / d ell_P = 0
d S       / d ell_P = -R_species ln(10)
d ell_S   / d ell_P = -R_species/S
dc1       / d ell_P = 0
dc2       / d ell_P = 0
```

Equivalently, in natural pressure:

```text
dS / d ln(P) = -R_species
dS / dP      = -R_species/P
```

The corresponding Python implementation is intentionally direct:

```text
parse_table(path)
  reads rows with header: logT, nps
  reads entries with:    logP, comp1, comp2, logrho, logS, logu

pressure_points(logp_min, first.logp)
  returns the new grid points from logp_min through first.logp - dlogP
  verifies that logp_min is on the dlogP = 0.20 SCVH grid

extend_row(row, logp_min)
  skips rows already extended to logp_min
  skips rows whose first point is still at logP ~= 4
  computes entropy0 = 10^first.logs
  computes rgas = 10^(first.logp - first.logrho - row.logt)
  estimates one-sided old-table slopes above first.logp
  chooses the C1 decay anchor from the ideal logRho=-15 point
  prepends ideal-reference entries with the C1 correction above the anchor
  adds an endpoint bump so the first one-cell slope is also continuous
  limits each composition pair to keep c1 + c2 <= 1

validate_rows(rows, logp_min)
  checks 63 logT rows
  checks dlogP spacing
  checks that extended rows start at logp_min

write_table(path, rows)
  writes the updated .asc.data file in the same 6-column format
```

### Confidence level for the low-density extension

Because the original first-stage script has not been found, the new
helper cannot be claimed to reproduce the historical preprocessing
exactly.  The defensible claim is narrower:

```text
below the existing processed edge, continue the already-packaged
low-pressure SCVH tail as a neutral/molecular ideal-gas boundary fill.
```

That is the physically expected limit as `P -> 0` at low temperature,
but it is still an extrapolated MESA preprocessing choice, not published
SCVH table coverage.

The helper is written to avoid changing historical processed values.  It
parses each existing entry with its raw input line and writes those lines
back unchanged.  Only the row header count and the new prepended
low-pressure entries change.

The local consistency checks for the prepended region are now:

```text
h_tab.asc.data:
  extended_rows = 25
  new_pressure_points = 675
  existing-entry diffs above logP=-0.60 = 0
  max one-step splice slope jump after discrete-slope correction:
     comp1 = 5e-6, comp2 = 5e-6, logu = 8.9e-15
  max current-minus-ideal at logRho=-15:
     comp1 = 2.88e-4, comp2 = 9.94e-5, logu = 2.01e-4

he_tab.asc.data:
  extended_rows = 25
  new_pressure_points = 675
  existing-entry diffs above logP=-0.60 = 0
  max one-step splice slope jump after discrete-slope correction:
     comp1 = 5.6e-16, comp2 = 1.4e-16, logu = 8.9e-15
  max current-minus-ideal at logRho=-15:
     comp1 = 3.21e-4, comp2 = 1.61e-4, logu = 3.67e-4
```

The one-cell slope through the splice is now continuous to printed-table
precision, while the density-floor values remain ideal to the same
precision.  These checks show that the new points follow the anchored
derivative-matched continuation, not that the missing historical
generator would have produced exactly the same numbers.

The main remaining risks are:

```text
1. one-sided derivatives at the old logP = -0.60 splice,
2. rows near logT ~= 4 where the processed edge is not a pure cold H2/He
   asymptote,
3. bicubic interpolation/smoothing behavior after the values are loaded,
4. the final eosDT OPAL/SCVH/HELM region selection after regeneration.
```

Those risks have to be checked in the regenerated `eosDT` table and the
plotter.  For the target cold low-density corner, the required checks
are that `mu`, `free_e`, and the neutral/molecular fractions follow the
neutral/molecular limit and that the selected EOS is the regenerated
OPAL/SCVH table rather than HELM fallback.

The old fraction-plot image artifacts were removed from `notes/SCVH_eos`.

This plot is made directly from:

```text
eos/eosDT_builder/eos_input_data/scvh/h_tab.asc.data
eos/eosDT_builder/eos_input_data/scvh/he_tab.asc.data
```

using the current `scvh_core.f90` interpretation:

```text
h_tab.asc.data:
  comp1 -> x_H2
  comp2 -> x_H
  x_H+  = 1 - (x_H2 + x_H)

he_tab.asc.data:
  comp1  -> x_He
  comp2  -> x_He+
  x_He++ = 1 - (x_He + x_He+)
```

With the derivative-matched helper, the new pressure entries are no
longer exactly pressure-flat right next to `logP = -0.60` in rows where
the old processed table has a nonzero composition slope.  They decay
back to the ideal frozen-composition tail toward `logP = -6.00`.
Thus the plotted variation is still mostly with `logT`, but the splice
neighborhood now carries the old one-sided pressure derivative.

The broader fraction diagnostic extended through the old processed tail
up to the original `logP = 4.00` edge.  It showed that the historical
`logP = -0.60` to `4.00` region is not pressure-flat: composition
fractions change with pressure, especially near `logT ~= 3.2` to `4.0`.
The local extension below `logP = -0.60` therefore should be understood
as a low-pressure boundary continuation from the existing processed edge,
not as a reconstruction of the full historical `logP = 4.00` to `-0.60`
processing.

An additional SVN archive search for the original extension equation did
not find a formula-bearing generator.  The fetched history contains:

```text
2007-01-29, SVN r119:
  imports h_tab.dat, h_tab_p1.dat, h_tab_p2.dat, he_tab.dat
  cold rows start at logP = 4.00

2010-04-15, SVN r2260:
  adds generated h_tab.asc.data and he_tab.asc.data
  cold rows start at logP = -0.60
  includes make_SCVH_plots.f and write_file_for_scvh
```

but no source routine was found that explicitly computes the
`logP = 4.00` to `-0.60` extension.  Thus the superseded table extension
should be described as inferred from the existing processed edge and the
ideal low-pressure limit, not copied from the original generator.

The superseded table-generation implementation extended 25
low-temperature rows in each file:

```text
logT = 2.10 through 4.02: nps = 126, logP = -6.00 ... 19.00
logT = 4.10 through 7.06: nps = 76,  logP =  4.00 ... 19.00
```

In that superseded attempt, the compressed builder input tarball was
repacked so the changed processed SCVH files were present in:

```text
eos/eosDT_builder/eos_input_data.tar.xz
```

The current branch has restored that tarball and the extracted SCVH
`.asc.data` files to the packaged state.

Superseded code changes made with the new processed grid:

```text
eos/eosDT_builder/src/scvh_core.f90
  NlogPs = 126
  logP_min = -6.00d0
  Rho_min = 1.0d-15
  write_file_for_scvh now computes nps from NlogPs/logP_max/dlogP

eos/eosDT_builder/src/create_eos_files.f90
  logRho_min = -15d0
  version_number = 52
  output rows use explicit whitespace between numeric fields

eos/eosDT_builder/src/eos_regions_defs.dek
  low-density OPAL/SCVH side moved to logRho ~= -15
  logQmin = -10.09d0

eos/defaults/eos.defaults
  logRho_min_OPAL_SCVH_limit = -15.0d0
  logQ_min_OPAL_SCVH = -10.09d0

eos/private/eosdt_eval.f90
  runtime low-density OPAL/SCVH tail limits moved below logRho = -15
```

No MESA compile or eosDT regeneration run was done after these source
and input-data edits.

## Self-consistency and derivative status of the superseded generated tail

The new `logP = -6.00` extension is not a recomputation of the published
SCVH EOS.  It is an ideal neutral/molecular pressure tail spliced onto
the already-processed MESA SCVH files.  The extension is value-level
self-consistent in the same limited sense as the historical low-pressure
tail:

```text
at fixed T:

d log10(rho) / d log10(P) = 1
u(P) = u(P0)
S(P) = S(P0) + R_species ln(P0/P)
```

where `P0` is the previous processed edge and `R_species` is inferred
from the existing edge value:

```text
R_species = P0/(rho0 T)
```

The implementation is in:

```text
eos/eosDT_builder/extend_scvh_low_p_tail.py
```

and the exact formulas are:

```text
logrho(logP) = logrho(logP0) + logP - logP0
u(logP)      = u(logP0)
S(logP)      = S(logP0) + R_species ln(P0/P)
```

The composition columns are held fixed below the previous edge.  In the
cold H and He rows, this preserves the already-processed neutral or
molecular state and avoids introducing a pressure-dependent ionization
closure in the artificial tail.

The derivatives should be described more carefully.  The current
`.asc.data` files contain only values:

```text
logP, comp1, comp2, logrho, logS, logu
```

They do not contain explicit analytic derivatives.  Current
`eos/eosDT_builder/src/scvh_core.f90` reads those values, calls
`fill_and_smooth`, and then obtains `d/dlogP` and `d/dlogT` from the
bicubic interpolation setup/evaluation:

```text
eos/eosDT_builder/src/scvh_core.f90
  read_file_for_scvh
  fill_and_smooth
  interp_mkbicub_db
  interp_evbicub_db
```

So the derivatives are locally derived from the extended value table and
are consistent with the interpolation machinery, but they are not an
independently verified analytic derivative table.  In the low-pressure
tail interior they should reflect the ideal-gas splice above.  Near the
old edge, near the high-temperature rows that are still filled/smoothed
in memory, and anywhere bicubic smoothing sees rounded table values, the
derivatives need to be validated by the builder/plotter output rather
than assumed exact.

The fetched 2010-04-15 SVN patch that first shows the processed
`logP = -0.60` `.asc.data` files does not include a preserved recipe or
comment explaining how the historical `logP = 4.00` to `-0.60` tail was
generated.  It does show that the active path switched to reading
`.asc.data`, filling/smoothing value arrays, and deriving partials from
the interpolation rather than reading pretabulated derivative columns.
Its comments are about smoothing and entropy of mixing; they are not a
validation note for the low-pressure-tail derivatives.

## Note on the higher-density plot artifact

The boxed feature near `logRho ~ -2.5` to `3` and `logT ~ 3` to `5` in the
mu plot is not in the low-pressure region extended by
`extend_scvh_low_p_tail.py`.  In pressure coordinates that region is
roughly the SCVH pressure-ionization part of the table.  For a neutral
ideal-gas estimate,

```text
logP ~= logRho + logT + log10(k_B / particle_mass)
```

so points around `logRho = 0`, `logT = 4` correspond to:

```text
H2: logP ~= 11.6
He: logP ~= 11.3
```

Across the full boxed range, a rough neutral-gas pressure estimate is:

```text
logRho = -2.5, logT = 3.0 -> logP ~= 7.8 to 8.1
logRho =  3.0, logT = 5.0 -> logP ~= 15.3 to 15.6
```

so the box crosses both the internal OPAL/SCVH transition and the
high-pressure SCVH ionization structure.  The historical
`extra_smoothing(compv1_he_si, ...)` call only targets the subset with
`logT < 4.5` and `logP > 10`.

That is close to the historical smoothing region called out in the
2010-04-15 SVN patch.  The r2260 patch added this active call:

```fortran
! do extra smoothing of compv1_he for logT < 4.5 and logP > 11.0 in
! this is where we get problems from the He+ to He++ pressure ionization
call extra_smoothing(compv1_he_si, 1, logTs(1), 4.5d0, 10d0, logPs(NlogPs))
```

Before this branch, the same call existed in the current MESA source but
was commented out:

```fortran
!call extra_smoothing(compv1_he_si,1,logTs(1),4.5d0,10d0,logPs(NlogPs))
```

Local old releases checked (`r11701` through `r15140`) also have that
call commented out.  A regenerated mu plot with this smoothing restored
did not remove the boxed feature.  The smoothing has therefore been
returned to the commented-out state:

```fortran
! Extra smoothing of the composition tables can move the ionization
! structure without making the thermodynamic table consistent.
!call extra_smoothing(compv1_he_si, 1, logTs(1), 4.5d0, 10d0, logPs(NlogPs))
```

The practical reason is that this filter smooths only `compv1_he`, not
the paired composition column or the thermodynamic quantities.  That can
move the ionization structure seen by `mu` without making `rho`, `S`,
`u`, and the composition table mutually consistent.  Since it did not
solve the regenerated plotter artifact, it should stay off unless a
separate, localized He pressure-ionization artifact is being isolated.

The low-pressure Python helper prepends points only from `logP = -6.00`
through the old `logP = -0.60` splice region.  It does not change the
`logP > 10` pressure-ionization values directly.  If the center-box
feature worsened after regeneration, the more likely causes are the
current `scvh_core` interpolation/smoothing state, the OPAL/SCVH/HELM
blend-region regeneration, or the builder fallback changes, not the
low-pressure continuation itself.

The current conclusion is that the cause is elsewhere in the regenerated
OPAL/SCVH table, the OPAL/SCVH/HELM region-selection path, or the
builder fallback path.  The plotter feature is a `rho,T`-plane transition
feature; it is not direct evidence that the `logP = -6.00` low-pressure
tail itself is wrong.

The regenerated plotter uses the runtime EOS path, not just the builder
region file.  The builder-side `eos_regions_defs.dek` now has the
low-density SCVH limits moved to the `rho = 1d-15` floor, but
`eos/private/eosdt_eval.f90` still contains the runtime OPAL/SCVH to
HELM transition logic with hard-coded values:

```fortran
logT6 = 4.890d0
logT5 = 4.899d0   ! problems with blend here so just jump
```

That runtime transition intersects the boxed `mu` feature in the
`rho,T` plot.  This is another reason the one-column SCVH composition
smoothing is the wrong knob for this artifact.

### Regenerated-table diagnostic for the boxed `mu` feature

The next diagnostic was done directly on the regenerated `eosDT` table
files, without compiling or running MESA.  For the plotter setup with
`Z = 0.02` and `X = 0.70`, the table comparison samples the two adjacent
composition tables as

```text
f(X = 0.70, Z = 0.02) ~= 0.5 f(X = 0.60, Z = 0.02)
                         + 0.5 f(X = 0.80, Z = 0.02)
```

and maps table coordinates with

```text
logQ = logRho - 2 logT + 12.
```

The old table-comparison image artifacts were removed from
`notes/SCVH_eos`.

The baseline was the local Git LFS object for the packaged
`eos/eosDT_data.tar.xz` table at `HEAD`, whose relevant files are version
51.  The regenerated files are version 53.

The boxed-region result is clear:

```text
current:  mesa-eosDT_02z60x.data, version 53
HEAD/LFS: mesa-eosDT_02z60x.data, version 51

current:  mesa-eosDT_02z80x.data, version 53
HEAD/LFS: mesa-eosDT_02z80x.data, version 51
```

Both versions have the same `logQ` and `logT` table shape:

```text
nQ = 527, logQ = -10.09 .. 5.69
nT = 306, logT =  2.10 .. 8.20
```

But the regenerated table has many abrupt jumps in the boxed region
where the packaged table is smooth in the pressure direction:

```text
X = 0.60, Z = 0.02, -3.5 <= logRho <= 4.5, 2.0 <= logT <= 5.0
  current dQ jumps with |Delta mu| > 0.4: 1854, max 1.79786
  HEAD    dQ jumps with |Delta mu| > 0.4:    0, max 0.15416

X = 0.80, Z = 0.02, same box
  current dQ jumps with |Delta mu| > 0.4: 1926, max 1.43022
  HEAD    dQ jumps with |Delta mu| > 0.4:    0, max 0.16259
```

The current table also writes a diagnostic `HELM` column from
`eos/eosDT_builder/src/create_eos_files.f90`.  That column is not read by
the runtime loader, which still has `num_eos_file_vals = 16` in
`eos/private/eosdt_load_tables.f90`, but it is useful for diagnosing how
the table was generated.  In the red-box region, the patchy `mu` and
`log_free_e` structure coincides with patchy builder `HELM_fraction`.

The expected smooth builder fraction is the region value `alpha` from
`eos/eosDT_builder/src/eos_regions_code.dek`:

```text
beta = 1 - alpha
```

where `alpha = 1` means HELM and `alpha = 0` means OPAL/SCVH.  In the
nominal pure-OPAL/SCVH part of the box, `alpha_region = 0`.  The
regenerated table nevertheless has many points where the written
`HELM_fraction` differs strongly from that smooth region value:

```text
X = 0.60, Z = 0.02
  points with |HELM_fraction - alpha_region| > 0.1: 2871
  points with |HELM_fraction - alpha_region| > 0.5: 1948
  points with current HELM_fraction = 1 and alpha_region < 0.1: 763

X = 0.80, Z = 0.02
  points with |HELM_fraction - alpha_region| > 0.1: 2942
  points with |HELM_fraction - alpha_region| > 0.5: 2028
  points with current HELM_fraction = 1 and alpha_region < 0.1: 827
```

This is the builder fallback path in
`eos/eosDT_builder/src/helm_opal_scvh_driver.f90`, not the runtime
plotter:

```fortran
call interpolate_opal_scvh(...)
call check_results
if (ierr /= 0) then
   ...
   if (ierr /= 0) then ! use helm instead
      beta = 0
      alfa = 1
   end if
end if
...
HELM_fraction = alfa
```

Thus the regenerated table can be in a nominal OPAL/SCVH region, fail
`check_results`, and be written as pure HELM at isolated points.  That
explains the speckled/patchy `mu` structure.  Smoothing only
`compv1_he_si` cannot fix this because the artifact is no longer just a
composition-column roughness; it is a discontinuous EOS-source fallback
in the generated table.

The existing SCVH input values above the old splice point were also
checked against the local Git LFS baseline for
`eos/eosDT_builder/eos_input_data.tar.xz`.  For both
`h_tab.asc.data` and `he_tab.asc.data`, all entries with
`logP >= -0.60` match the packaged `HEAD` input exactly.  The local
extension added 675 lower-pressure points per species and did not alter
the original entries:

```text
h_tab.asc.data:  new points below -0.60 = 675, existing-entry diffs = 0
he_tab.asc.data: new points below -0.60 = 675, existing-entry diffs = 0
```

The noisy regenerated `mesa-eosDT_02z60x.data` and
`mesa-eosDT_02z80x.data` files being inspected here were written around
13:32-13:33 on 2026-05-13.  The derivative-matched helper and extracted
SCVH inputs were updated later, around 14:48 on 2026-05-13.  Therefore
the noisy plot is evidence against the earlier value-only low-pressure
prepend, not against the current derivative-matched files.

The old value-only tarball had large first-derivative jumps at the
`logP = -0.60` splice because the new tail was flat in composition and
`logu` while the old processed table above the splice was not:

```text
old tarball, one-step slope jump across logP = -0.60
  H  comp1  max |above - below| = 0.215888
  H  comp2  max |above - below| = 0.193234
  H  logu   max |above - below| = 0.2555
  He comp1  max |above - below| = 0.347985
  He comp2  max |above - below| = 0.173995
  He logu   max |above - below| = 0.388
```

The current extracted files, generated by
`eos/eosDT_builder/extend_scvh_low_p_tail.py`, match the one-cell
pressure slope across the splice while forcing the correction to decay
by the ideal `rho = 1d-15` anchor:

```text
current extracted files, one-step slope jump across logP = -0.60
  H  comp1  max |above - below| = 5e-6
  H  comp2  max |above - below| = 5e-6
  H  logrho max |above - below| = 8.9e-15
  H  logs   max |above - below| = 8.9e-15
  H  logu   max |above - below| = 8.9e-15
  He comp1  max |above - below| = 5.6e-16
  He comp2  max |above - below| = 1.4e-16
  He logrho max |above - below| = 8.9e-15
  He logs   max |above - below| = 8.9e-15
  He logu   max |above - below| = 8.9e-15
```

At the density floor itself, the current files reproduce the
ideal-reference tail to table precision:

```text
current extracted files, interpolated current-minus-ideal at logRho=-15
  H  comp1  max |Delta| = 2.88e-4
  H  comp2  max |Delta| = 9.94e-5
  H  logrho max |Delta| = 1.10e-4
  H  logs   max |Delta| = 7.66e-5
  H  logu   max |Delta| = 2.01e-4
  He comp1  max |Delta| = 3.21e-4
  He comp2  max |Delta| = 1.61e-4
  He logrho max |Delta| = 9.18e-5
  He logs   max |Delta| = 5.32e-5
  He logu   max |Delta| = 3.67e-4
```

The old slope-diagnostic image artifact was removed from `notes/SCVH_eos`.

The tarball `eos/eosDT_builder/eos_input_data.tar.xz` was then rebuilt
from the current extracted `eos_input_data` directory at 15:15 on
2026-05-13 and verified byte-identical for:

```text
eos_input_data/scvh/h_tab.asc.data
eos_input_data/scvh/he_tab.asc.data
```

The next rebuild should therefore test the derivative-matched SCVH input.
If the speckled `mu` feature remains after that rebuild, the culprit is
not the old value-only splice.  The next place to inspect is then the
builder fallback path in `helm_opal_scvh_driver.f90`: isolated
`check_results` failures in a nominal OPAL/SCVH region still write pure
HELM points and will appear as patchy `mu`/`free_e`.

Rebuild checks:

```text
1. Regenerate the builder tables from the synced input tar/directory.
2. Confirm mesa-eosDT_02z60x.data and mesa-eosDT_02z80x.data have mtimes
   later than the 15:15 input-tar sync.
3. Refresh the runtime eosDT data tar before using the plotter, otherwise
   the plotter may still read stale 13:32 tables.
4. Replot mu and compare the diagnostic HELM_fraction with the smooth
   alpha_region from eos_regions_code.dek in the red-box region.
```

### Fix direction

If derivative matching is not enough, the fix should make the
low-pressure extension local.  The old
`logP >= -0.60` SCVH surface should be evaluated with the same effective
grid and smoothing/fill rules as the packaged version, while the new
`logP < -0.60` tail should be handled by a separate branch.

The important boundary is:

```text
logP_edge = -0.60
```

For the old domain:

```text
logP >= logP_edge:
   use the original 99-point SCVH pressure grid
   logP = -0.60, -0.40, ..., 19.00
```

For the new tail:

```text
logP < logP_edge:
   use only the low-pressure continuation
```

The tail equations remain the edge-matched ideal continuation:

```text
R_species = P0 / (rho0 T)
logrho(P) = logrho0 + logP - logP0
u(P) = u0
S(P) = S0 + R_species ln(P0/P)
composition(P) = composition0
```

This avoids letting the added columns at `logP = -6.00 .. -0.80`
participate in the bicubic spline that is used for the historical
pressure-ionization region.  In code terms, the least invasive clean
approach is to split the SCVH interpolation:

```text
if logP >= -0.60:
   evaluate a 99-column "old" bicubic table
else:
   evaluate the low-pressure continuation branch
```

An intermediate diagnostic guard is also useful while regenerating:

```text
alpha_region = value from eos_regions_code.dek
HELM_fraction = alfa written by helm_opal_scvh_driver.f90
```

If `abs(HELM_fraction - alpha_region)` is large in a nominal pure
OPAL/SCVH region, the builder should report it instead of silently
writing scattered HELM fallback values.  This would have caught the
boxed feature before the plotter stage.

## Consistency with the other branch changes

The branch has several changes that are not part of the raw pressure-tail
generation but are consistent with making the extended tail usable.

### Builder source names

```text
eos/eosDT_builder/Makefile
```

The builder `Makefile` now points at the `.f90` source names:

```text
src/azbar.f90
src/create_eos_files.f90
src/helm_opal_scvh_driver.f90
src/opal_scvh_driver.f90
src/scvh_core.f90
```

This matches the current source tree.  It is a build-system consistency
fix, not an EOS-physics change.

The generated private-source rules also use directory prerequisites of
the form:

```make
$(BUILD_DIR_MODULE)/private/.
```

This matches the top-level MESA `make/setup-builddir.mk` directory
target pattern.  Without the trailing `/.`, `make` can fail with a
message like:

```text
No rule to make target '../../build/eos-dt-builder/private',
needed by '../../build/eos-dt-builder/private/helm_alloc.f90'
```

### eosDT mode rather than standalone SCVH mode

```text
eos/eosDT_builder/src/create_eos_files.f90
```

The branch keeps:

```fortran
scvh_only = .false.
```

This is needed for the normal packaged `eosDT` product.  With
`scvh_only = .true.`, the builder writes standalone `eosSCVH_data`
instead of `eosDT_data`, which is not the table family used by the MESA
default OPAL/SCVH/eosDT path.

The data-row output format was also widened and given explicit
whitespace separators.  The previous fixed-width format could produce
glued values when a dimensionless result became large, for example:

```text
-0.762701139.11869
```

which the runtime loader cannot parse as separate `log_free_e` and
`gamma1` values.  This is a table-formatting issue in the builder output,
not by itself evidence that the SCVH extension failed.

### Builder-side HELM avoidance

```text
eos/eosDT_builder/src/helm_opal_scvh_driver.f90
```

The builder now computes:

```fortran
logdin = logRho + log10(max(1d-99, zbar/abar))
logQ = logRho - 2*logT + 12d0
avoid_helm = (logT <= logT_HELM .or. logdin <= eos_ht% logdlo .or. logQ <= logQmin)
```

and forces `alfa = 0` when `avoid_helm` is true.  This makes the builder
prefer the OPAL/SCVH side in the cold low-density region instead of
falling through to HELM, which is the branch behavior wanted for the
newly extended SCVH tail.

The new `get_ideal_tail` in the same file is a safety fallback, not the
primary fix.  It is used when HELM should be avoided and OPAL/SCVH still
fails.  It prevents a hard HELM failure or a fully ionized HELM value
from being written by the builder in that corner.  It does not replace
the SCVH low-pressure tail and it is not a molecular-H2 reconstruction;
the molecular/neutral cold behavior is expected to come from the
extended processed SCVH input.

### Builder region limits

```text
eos/eosDT_builder/src/eos_regions_defs.dek
```

The builder-side OPAL/SCVH-to-HELM region constants were moved from the
old low-density edge to the new one:

```text
logRho5 = -15.0
logRho6 = -15.1
logRho7 = -15.5
logRho8 = -15.9
logQmin = -10.09
```

These values are consistent with `scvh_core` using `Rho_min = 1d-15`
and with `create_eos_files` sampling down to `logRho_min = -15`.

### Runtime OPAL/SCVH selection limits

```text
eos/defaults/eos.defaults
eos/private/eosdt_eval.f90
```

The runtime OPAL/SCVH limits were also moved:

```text
logRho_min_OPAL_SCVH_limit = -15.0
logQ_min_OPAL_SCVH = -10.09
```

and the hard-coded low-density tail markers in `eosdt_eval.f90` now sit
below `logRho = -15`:

```text
logRho5 = -15.0
logRho6 = -15.001
logRho7 = -15.90
logRho8 = -15.99
```

This is needed so the regenerated OPAL/SCVH/eosDT table can actually be
selected at the new cold low-density edge.  For example:

```text
logT = 3.00
logRho = -15.00
logQ = logRho - 2 logT + 12 = -9.00
```

which is now above the runtime `logQ_min_OPAL_SCVH = -10.09` cutoff.

The runtime HELM/ideal fallback in `get_level6_for_eosdt` has not been
changed in this branch.  The intended behavior at the target point is
therefore to avoid reaching level 6 by selecting the regenerated
OPAL/SCVH table.  If OPAL/SCVH is disabled, outside composition limits,
or not regenerated/installed, the old HELM/ideal fallback policy can
still be reached.

### Plotter changes

```text
eos/plotter/inlist_plotter
eos/plotter/src/eos_plotter.f90
```

The plotter inlist now selects `i_var = 4`, which is `i_mu`, instead of
the special region-map mode `i_var = 0`.  That is a diagnostic choice for
checking the molecular/neutral versus ionized `mu` behavior.

The plotter source change keeps non-blend regions masked when
`only_blend_regions` is enabled without skipping the rest of the point
logic.  This is a plotting/diagnostic change only; it does not affect the
EOS tables or runtime EOS selection.

## Superseded remaining required data step

The source and builder-input changes are not enough by themselves for a
normal MESA run to see the new values.  The final packaged table still
has to be regenerated and exported:

```text
eos/eosDT_data.tar.xz
```

and then installed into the runtime data directory.  Until that is done,
the plotter or MESA may still read the old packaged `eosDT` data and
show the old HELM-like cold-corner behavior.

## Verification that the low-pressure tail is ideal gas

For molecular hydrogen, the ideal-gas relation is:

```text
P = rho k_B T / (2 m_H)
```

or:

```text
logrho = logP - logT - log10(k_B/(2 m_H))
```

Using cgs constants:

```text
log10(k_B/(2 m_H)) ~= 7.6154
```

At the historical pre-extension processed-table low-pressure edge:

```text
logT = 2.10
logP = -0.60

logrho(H2) = -0.60 - 2.10 - 7.6154
            = -10.3154
```

That matches the processed hydrogen low-pressure tail.

For neutral helium:

```text
P = rho k_B T / m_He
```

or:

```text
logrho = logP - logT - log10(k_B/m_He)
```

Using cgs constants:

```text
log10(k_B/m_He) ~= 7.3175
```

At the same historical pre-extension processed-table edge:

```text
logT = 2.10
logP = -0.60

logrho(He) = -0.60 - 2.10 - 7.3175
            = -10.0175
```

That matches the processed helium low-pressure tail.

At the new local edge:

```text
logT = 2.10
logP = -6.00

logrho(H2) = -6.00 - 2.10 - 7.6154
            = -15.7154

logrho(He) = -6.00 - 2.10 - 7.3175
            = -15.4175
```

Those values match the first rows in the repacked processed SCVH input
tarball.

The internal energies also match the same ideal-gas interpretation:

```text
H2: u ~= (5/2) k_B T / (2 m_H)
He: u ~= (3/2) k_B T / m_He
```

So below the original `logP = 4` edge, the MESA low-pressure tail should
be described as a neutral/molecular ideal-gas continuation attached to
the SCVH table, not as original SCVH table coverage.

## Why `rho = 1d-10` is not an original SCVH rho boundary

The original SCVH tables are not bounded by a single value of `rho`.
They are bounded in `logP` and `logT`, and the implied density varies
with temperature.

For example, at the original low-pressure edge `logP = 4`:

```text
H,  logT = 2.10 -> logrho = -5.7154
He, logT = 2.10 -> logrho = -5.4175
```

but at the hot end:

```text
H,  logT = 7.06 -> logrho = -11.2782
He, logT = 7.06 -> logrho = -10.8546
```

Thus the original tables do contain densities below `1d-10` in the hot,
low-pressure part.  The pre-extension `rho = 1d-10` value in MESA was a
runtime guard in the SCVH interface, not an original SCVH table floor.

## Extending the MESA tail toward `rho = 1d-15`

It is technically coherent to extend the existing MESA approach:

```text
old processed tail: logP = 4.0 down to -0.6
new processed tail: logP = 4.0 down to about -6.0
old runtime guard:  rho = 1d-10
new runtime guard:  rho = 1d-15
```

This would not recover additional published SCVH coverage.  It would be
a deeper MESA ideal-gas low-pressure continuation.

The pressure required for `rho = 1d-15` at the coldest SCVH grid point
is approximately:

```text
logP = logrho + logT + log10(k_B / particle_mass)
```

For neutral He at `logT = 2.10`:

```text
logP ~= -15 + 2.10 + 7.3175
     ~= -5.58
```

For molecular H2 at `logT = 2.10`:

```text
logP ~= -15 + 2.10 + 7.6154
     ~= -5.28
```

So `logP_min = -6.0` gives margin for the coldest low-density case.

With the existing `dlogP = 0.20` and `logP_max = 19`, the SCVH pressure
grid would need:

```text
NlogPs = (19 - (-6))/0.20 + 1
       = 126
```

## Superseded code and data changes implied by `logP_min = -6`

A consistent implementation would need all of the following:

1. Extend the processed SCVH data files.

   ```text
   eos/eosDT_builder/eos_input_data/scvh/h_tab.asc.data
   eos/eosDT_builder/eos_input_data/scvh/he_tab.asc.data
   ```

   For the cold rows currently carrying the full low-pressure extension,
   the number of pressure points would increase from `99` to `126`.
   High-temperature rows that currently start at `logP = 4` could still
   start at `logP = 4`, with `k_for_MIN_logP` marking the missing lower
   pressure region for the fill/smooth logic.

   The first-stage source files, if kept, would also need to be extended:

   ```text
   eos/eosDT_builder/eos_input_data/scvh/h_tab.dat
   eos/eosDT_builder/eos_input_data/scvh/he_tab.dat
   ```

   The current runtime does not read these files, but they are the
   closest available record of the historical preprocessing before the
   `.asc.data` files were produced.

2. Update the hard-coded SCVH pressure grid in `scvh_core.f90`.

   ```fortran
   integer, parameter :: NlogPs = 126, NlogTs = 63
   double precision, parameter :: dlogP = 0.20d0
   double precision, parameter :: logP_min = -6.00d0
   double precision, parameter :: logP_max = 19d0
   ```

3. Lower the SCVH runtime density guard.

   ```fortran
   parameter(Rho_min=1.0d-15)
   ```

4. Update any writer/helper assumptions in `scvh_core.f90`.

   The disabled writer previously used `nps = 99` for low-temperature
   rows and `nps = 76` for high-temperature rows.  If it is used to
   regenerate `.asc.data`, the low-temperature count needs to follow
   `NlogPs`.

5. Ensure the final `eosDT` builder samples the target low densities.

   `create_eos_files.f90` has a separate `logRho_min` floor.  If the
   packaged `eosDT` tables should include `logrho = -15`, that builder
   floor needs to be lowered consistently.

6. Ensure the OPAL/SCVH versus HELM/ideal region logic uses the intended
   EOS in the cold low-density region.

   Extending the SCVH source data alone is not sufficient if the final
   packaged `eosDT` region logic immediately selects HELM or another
   fallback there.

7. Regenerate the packaged `eosDT` data.

   The shipped `eosDT_data.tar.xz` contains precomputed values.  Runtime
   code changes alone do not alter those packaged tables.

## Wording to use

Good:

```text
MESA's processed SCVH files contain an ideal neutral/molecular
low-pressure extension below the published SCVH `logP = 4` edge,
historically down to `logP = -0.60`; this branch extends the processed
tail to `logP = -6.00`.
```

Good:

```text
The historical `rho = 1d-10` clamp is a MESA runtime guard in the SCVH
interface, not an original SCVH table boundary.
```

Good:

```text
The historical `rho = 1d-10` guard and the `logP = -0.60`
processed-table edge are numerically consistent, but the repository does
not prove which choice motivated the other.
```

Avoid:

```text
SCVH has an original density floor at rho = 1d-10.
```

Avoid:

```text
MESA extends SCVH to lower density by using more SCVH table data.
```

The lower-pressure region is better described as a MESA ideal-gas
continuation attached to the SCVH input.

## 2026-05-13 builder outcome

The current implementation no longer rewrites
`eos/eosDT_builder/eos_input_data/scvh/h_tab.asc.data` or
`eos/eosDT_builder/eos_input_data/scvh/he_tab.asc.data`.  Those input
tables are kept at the processed `logP_min = -0.60` edge.  The low
pressure extension is applied inside
`eos/eosDT_builder/src/scvh_core.f90` when the SCVH density solve asks for
`logP < logP_min`.

For each species, the continuation is anchored to the restored table edge
at fixed `logT`:

```text
logrho_ideal = logrho0 + logP - logP0
logu_ideal   = logu0
Rgas         = 10^(logP0 - logrho0 - logT)
S_ideal      = S0 + Rgas ln(10) (logP0 - logP)
```

The code blends from the edge derivative to the ideal derivative over one
pressure interval using a C1 correction

```text
x = (logP - logP_decay)/(logP0 - logP_decay)
g = x^3 - x^2
value = ideal_value + M g
M = (dvalue/dlogP|edge - dvalue/dlogP|ideal,edge) (logP0 - logP_decay)
```

so the interpolated value and the `logP` derivative are continuous at the
processed SCVH edge.  Composition fractions are edge-matched and limited
to remain in `[0,1]`.

Two builder-side issues appeared during regeneration:

1.  The internal OPAL/SCVH blend sometimes samples OPAL at very low
    density where OPAL returns `info=0` but derivative quantities are not
    finite.  `eos/eosDT_builder/src/opal_scvh_driver.f90` now treats
    nonfinite OPAL results like an OPAL miss when SCVH is already part of
    the blend, and keeps the SCVH value rather than falling through to
    HELM.

2.  The rectangular eosDT support table includes high-temperature,
    low-`Q` cells that runtime selection does not intend OPAL/SCVH to own.
    `eos/eosDT_builder/src/helm_opal_scvh_driver.f90` now retries HELM
    without electron-positron terms if the full HELM call fails.  If HELM
    still cannot fill that support cell, it uses a simple ideal gas plus
    radiation value for the HELM-side fill.  Runtime defaults in
    `eos/defaults/eos.defaults` and `eos/private/eosdt_eval.f90` were not
    changed for this.

The generated text format in `create_eos_files.f90` now writes
`gamma1`, `gamma3`, and `grad_ad` in exponential format so large support
values do not produce `************` fields that the table reader cannot
parse.

Validation run:

```text
cd eos/eosDT_builder
make
cp ../../build/eos-dt-builder/bin/eos-dt-builder ceos
MESA_DIR=$PWD/../.. ./rn
MESA_DIR=$PWD/../.. ./build_eosDT_4_export
cp data/eosDT_data/*.data ../../data/eosDT_data/

cd ../plotter
make clean
MESA_DIR=$PWD/../.. make plot
MPLBACKEND=Agg ./regions.py
MPLBACKEND=Agg ./plotter.py
```

`./rn` completed for `Z = 0.00`, `0.02`, and `0.04`; the exported tarball
is `eos/eosDT_data.tar.xz`.  The `make plot` target reaches the Python
plotting step, but `plotter.py` calls `plt.show()`, so noninteractive
verification should run the Python scripts with `MPLBACKEND=Agg`.

## 2026-05-13 remaining mu artifact

The regenerated `mu` plot still showed narrow noisy structure in the
OPAL/SCVH-owned region.  A region diagnostic showed that the
low-density, low-temperature feature is inside the OPAL/SCVH region
rather than a FreeEOS or runtime HELM boundary.

The support tables already contain isolated HELM fallbacks at the
affected points.  For example, near the processed SCVH pressure edge,

```text
logRho = -11.45, logT = 3.30, logQ = -6.05

mesa-eosDT_02z60x.data: mu = 0.668, HELM = 1
mesa-eosDT_02z80x.data: mu = 2.292, HELM = 0
```

and on the right side of the OPAL/SCVH island,

```text
logRho = 0.00, logT = 4.00, logQ = 4.00

mesa-eosDT_02z60x.data: mu = 0.668, HELM = 1
mesa-eosDT_02z80x.data: mu = 1.383, HELM = 0
```

Thus the artifact is not caused simply by runtime interpolation choosing
the wrong table cells.  The discontinuity is baked into the generated
eosDT support files before the runtime interpolant sees them.  Runtime
interpolation then blends across X, logQ, and logT and makes the isolated
fallback points visible as noisy structure.

The fallback was triggered by `check_results` in
`eos/eosDT_builder/src/helm_opal_scvh_driver.f90`, which rejected finite
OPAL/SCVH support values when `gamma1 >= 2` or `grad_ad >= 50`.  Those
large derivatives can occur near pressure-ionization boundaries.  The
builder now rejects nonfinite or nonpositive OPAL/SCVH thermodynamics,
but no longer converts otherwise finite support points to HELM only
because `gamma1` or `grad_ad` is large.  This keeps the eosDT support
table on the same EOS across the boundary instead of planting isolated
HELM points in the OPAL/SCVH region.

After rebuilding with that change, the right-side sample points at
`logRho = -0.57, logT = 4.18` and `logRho = 0.00, logT = 4.00` no longer
had `HELM = 1`, but the plotted `mu` features remained.  A second small
change in `eos/eosDT_builder/src/opal_scvh_driver.f90` broadened the
existing `search_for_SCVH` retry so all negative SCVH rejections try a
slightly smaller `logQ` before falling back.  This removed only a few
low-temperature fallback cells.

Remaining table diagnostic after the second rebuild:

```text
left window: -13 < logRho < -8, 2.8 < logT < 3.6

mesa-eosDT_02z60x.data: 47 cells with HELM > 0.5
mesa-eosDT_02z80x.data: 55 cells with HELM > 0.5

example:
logRho = -11.10, logT = 3.16, logQ = -5.42
mesa-eosDT_02z60x.data: mu = 2.200, HELM = 0
mesa-eosDT_02z80x.data: mu = 0.573, HELM = 1
```

The current interpretation is that the plotted features are still caused
by eosDT support cells that fall back to HELM inside or next to the
OPAL/SCVH-owned support region.  The fallback is now mostly from SCVH
evaluation failures, not from the finite-derivative guard in
`check_results`.  The next clean fix should avoid writing fully-ionized
HELM values into low-temperature OPAL/SCVH support cells; either the SCVH
evaluation needs a more robust local fill there, or the builder needs a
neutral/SCVH-like analytic support fill rather than a HELM fallback.

## 2026-05-13 OPAL/SCVH-only support rebuild

The right-side feature was traced to using the full outer
HELM/OPAL-SCVH blend while generating `mesa-eosDT_*.data`.  Setting
`scvh_only = .true.` was not the right correction because it changes the
output path to `data/eosSCVH_data/mesa-SCVH_*.data`.  The cleaner
builder setting for the eosDT support files is

```fortran
opal_scvh_only = .true.
scvh_only = .false.
```

in `eos/eosDT_builder/src/create_eos_files.f90`.  This keeps the output
in `data/eosDT_data` while preventing the outer HELM blend from being
baked into the OPAL/SCVH support table.

After rebuilding with this setting, right-side support cells no longer
carry HELM fraction:

```text
right window: -2 < logRho < 2, 3.4 < logT < 4.8

mesa-eosDT_02z60x.data: 0 cells with HELM > 0.1
mesa-eosDT_02z80x.data: 0 cells with HELM > 0.1

example:
logRho = 1.00, logT = 3.80, logQ = 5.40
mesa-eosDT_02z60x.data: mu = 0.835, HELM = 0
mesa-eosDT_02z80x.data: mu = 0.668, HELM = 0
```

Any remaining right-edge shape in the blended runtime plot is therefore
from the runtime EOS-region boundary and interpolation through the
OPAL/SCVH support values, not from HELM-filled support cells in the
regenerated table.

The saved pre-change reference plot was compared against the current
`opal_scvh_only` rebuild.  The right side is no longer affected by
generated HELM fractions, but the current rebuild still does not exactly
reproduce the old right-edge appearance.  That remaining change is
therefore tied to accepting different OPAL/SCVH support values near the
pressure-ionization/high-`Q` edge, not to the analytic low-`P` tail.

## 2026-05-13 FreeEOS low-temperature fallback row

The horizontal feature in the plotter's `lnE` partial with respect to
`lnT` at low density and `logT = 3.0..3.1` is not from the SCVH analytic
tail.  It coincides with the runtime FreeEOS lower-temperature blend:

```text
eos/defaults/eos.defaults
logT_min_FreeEOS_lo = 3.0d0
logT_min_FreeEOS_hi = 3.1d0
```

The FreeEOS support tables also contain fallback rows at the same edge.
For example, near `logRho = -15`:

```text
data/eosFreeEOS_data/mesa-FreeEOS_02z70x.data
logT = 3.02, MESA = 0, logE column = 12.9630
logT = 3.00, MESA = 1, logE column = 27.4192
logT = 2.98, MESA = 1, logE column = 27.2446
```

The `MESA = 1` rows come from
`eos/eosFreeEOS_builder/src/free_eos_table.f90:mesa_eos_eval`.  That
routine copies `eosDT_get` results into the FreeEOS table fallback row.
The first three logarithmic table columns are written as base-10 logs by
`free_eos_eval`, but `eosDT_get` returns `lnPgas`, `lnE`, and `lnS`.
The fallback rows therefore put natural-log values into columns that
`eos/private/eosdt_load_tables.f90` later reads as base-10 and multiplies
by `ln10`.

That unit mismatch makes the fallback row discontinuous with the first
valid FreeEOS rows.  The runtime blend derivative then includes the
expected term

```text
d alpha/d lnT * (lnE_FreeEOS - lnE_lower_priority)
```

which turns the fallback discontinuity into a horizontal stripe in the
`lnE` partial.  The clean fix is to correct the fallback row units in the
FreeEOS builder and regenerate the FreeEOS data.  The code change is in
`eos/eosFreeEOS_builder/src/free_eos_table.f90:mesa_eos_eval`: after
copying `eosDT_get` results into the FreeEOS-format row, convert
`i_lnPgas`, `i_lnE`, and `i_lnS` from natural log to base-10 log before
writing.  The FreeEOS table version was also bumped from 51 to 52 in the
builder defaults and shipped inlists so regenerated tables invalidate old
binary caches.  This is separate from the OPAL/SCVH low-`P` tail and from
the generated OPAL/SCVH eosDT support tables.

The literal FreeEOS builder rebuild was blocked in this checkout because
the external `free-eos` package was not available through `pkg-config`.
The existing extracted FreeEOS tables were therefore corrected
mechanically in the same way the fixed builder writes them: for rows with
`MESA = 1`, divide only the `logPgas`, `logE`, and `logS` table columns
by `ln10`, then change the table version header from 51 to 52.  The
runtime data in `data/eosFreeEOS_data` and the shipped archive
`eos/eosFreeEOS_data.tar.xz` were updated from those corrected tables.

For the earlier example near `logRho = -15`, the corrected fallback rows
are now

```text
data/eosFreeEOS_data/mesa-FreeEOS_02z70x.data
logT = 3.02, MESA = 0, logE column = 12.9332
logT = 3.00, MESA = 1, logE column = 11.9080
logT = 2.98, MESA = 1, logE column = 11.8322
```

The regenerated plotter image for `i_lnE`, partial with respect to
`lnT`, no longer shows the severe black horizontal fallback-row stripe at
`logT = 3.0`.  The diagnostic images are saved as:

```text
notes/SCVH_eos/eos_plotter_lnE_dlnT_v52_freeeos_fix.png
notes/SCVH_eos/eos_plotter_mu_current_v52_freeeos_fix.png
notes/SCVH_eos/eos_plotter_frac_HELM_current.png
notes/SCVH_eos/eos_plotter_frac_OPAL_SCVH_current.png
```

Those fraction plots show a separate issue: the visible feature on the
right side of the SCVH island is the runtime OPAL/SCVH high-`Q` boundary
falling back to HELM (`i_frac_OPAL_SCVH = 0`, `i_frac_HELM = 1`).  That
is controlled by `eos/private/eosdt_eval.f90:get_opal_scvh_alfa_and_partials`
and the OPAL/SCVH controls in `eos/defaults/eos.defaults`, not by the
FreeEOS fallback-row unit bug.
