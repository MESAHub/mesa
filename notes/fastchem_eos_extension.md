# FastChem EOS extension for MESA

Status: note only.  No MESA compile, MESA run, FreeEOS table generation, or
FastChem table generation has been performed for this note.

## Question

Can MESA be extended with a FastChem-based cool-chemistry EOS table, using a
free-energy/thermochemistry construction for self-consistent thermodynamics,
while also understanding how this overlaps with regenerated FreeEOS and the
current OPAL/SCVH-like regime?

## Short answer

FreeEOS can probably be regenerated over the OPAL-like regime, and that is
the least controversial extension target.  MESA already describes its shipped
FreeEOS tables as covering a domain similar to OPAL, and the current MESA
FreeEOS table grid is broad:

```text
2.1 <= log10(T) <= 8.2
-10.09 <= log10(Q) <= 5.69
Q = rho / T^2 * 1d12
```

Equivalently,

```text
log10(rho) = log10(Q) + 2*log10(T) - 12
```

so the full grid envelope is approximately:

```text
-17.9 <= log10(rho/[g cm^-3]) <= 10.1
```

This envelope is not a rectangular region in `logT,logRho`; the active limit
at fixed temperature is set by the `logQ` lines.  For example:

```text
log10(T) = 6
log10(rho)_min = -10.09 + 2*6 - 12 = -10.09
log10(rho)_max =   5.69 + 2*6 - 12 =   5.69
```

Going farther than the OPAL-like regime is a physics/validation problem, not
just a table-generation problem.  The high-density, low-temperature pressure
ionization region is exactly where Irwin reported multivalued FreeEOS
solutions and plasma phase transition discontinuities in the EOS1
free-energy model.

## Upstream FreeEOS status

MESA documentation and the current `eosFreeEOS_builder` path describe the
MESA tables as based on FreeEOS 2.2.1 with the EOS1 option suite.

Upstream FreeEOS has been updated since then.  The current SourceForge file
listing gives `free_eos-3.0.0.tar.gz` as the latest released source package,
modified 2022-06-05.  SourceForge project metadata also shows a later project
`Last Update` date, but the released source package list still shows 3.0.0 as
the latest release.

Important 3.0.0 changes:

- converted source from fixed-form Fortran 77 style to free-form `.f90`
  Fortran 2008;
- fixed many bugs and improved reliability/convergence behavior;
- updated physical constants, ionization potentials, and isotopic masses;
- added modern CMake support and documentation;
- reduced serious floating-point exceptions in the upstream tests.

Important limitation: the 3.0.0 release notes say the free-energy
semi-empirical models were not changed relative to 2.2.1.  Therefore 3.0.0 is
the right starting point for a cleaner generator, but it should not be assumed
to fix the EOS1 pressure-ionization model deficiencies.

## Free-energy tabulation idea

For fixed composition and fixed independent variables `(rho,T)`, FreeEOS is
already solving a free-energy minimization problem:

```text
minimize F(T,V,{N_i})
```

subject to abundance conservation and charge neutrality.  Once a consistent
minimum is found, thermodynamic quantities are derivatives of one potential:

```text
P = - (dF/dV)_{T,N}
S = - (dF/dT)_{V,N}
E = F + T*S
```

A table based on a single smooth interpolant for the specific Helmholtz free
energy,

```text
f = F / m
```

would be thermodynamically cleaner than independently interpolating many
derived quantities.  In a `logT,logRho` table one would need a consistent
mapping from derivatives in the interpolation variables to derivatives in
`T,V`, for example:

```text
rho = m/V
d/dV|T = -rho/V * d/drho|T
P = rho^2 * d f / d rho |_T
S = - d f / d T |_rho
E = f + T*S
```

This design could reduce internal inconsistency from interpolation.  It would
not remove discontinuities or branch ambiguity in the underlying FreeEOS
solution.

The current MESA FreeEOS builder is not organized this way.  It calls
FreeEOS, receives `logf`, and then tabulates the returned derived quantities:

```text
P, E, S, chiRho, chiT, Cp, Cv, dE/dRho, dS/dT, dS/dRho, mu, ...
```

The `logf` value appears in the current call path but is not written to the
MESA FreeEOS tables:

```text
eos/eosFreeEOS_builder/src/free_eos_table.f90
```

Specifically, `free_eos_eval` declares `logf`, passes it to `free_eos`, and
then fills the `result` array from the derived outputs.  The later
`write_table` path writes the derived outputs, not `logf`.

So a free-energy-based table would be a real structural change:

1. store `f` or `logf` as the primitive table quantity;
2. interpolate that primitive with a differentiable interpolation scheme;
3. derive `P`, `S`, `E`, `Cv`, `chiRho`, `chiT`, `Gamma1`, and related
   quantities from derivatives of that same interpolant;
4. use directly returned FreeEOS derivatives as validation data, not as
   independently interpolated primitives.

This would fix a specific class of EOS problems: interpolation-induced
thermodynamic inconsistency, for example where independently interpolated
`P`, `E`, and `S` no longer satisfy the Maxwell identities or the same
underlying free-energy differential.

It would not fix:

```text
FreeEOS convergence failure at a point
missing low-temperature species or chemistry
wrong pressure-ionization model physics
multivalued branches in the free-energy model
first-order discontinuities / plasma phase transitions
bad behavior from choosing the wrong branch
```

Thus a free-energy table can make the tabulated EOS more internally
consistent, but it cannot make an invalid or incomplete FreeEOS solution
physically correct.

## Pressure-ionization risk

Irwin's convergence paper reports that the more robust FreeEOS convergence
method reached conditions relevant to models between `0.1 Msun` and the
hydrogen-burning limit near `0.07 Msun`, plus hot brown-dwarf models just
below that limit.  However, the same paper reports EOS discontinuities in the
high-density, low-temperature pressure-ionization region.

The reported critical points for the H-related plasma phase transition are
roughly:

```text
log10(T/[K]) ~= 4.7 to 4.8
log10(rho_SI) ~= 3.0 to 3.2
```

The paper uses SI density units, so:

```text
rho_SI [kg m^-3] = 1000 * rho_cgs [g cm^-3]
log10(rho_cgs) = log10(rho_SI) - 3
```

Thus the transition neighborhood is approximately:

```text
log10(rho/[g cm^-3]) ~= 0.0 to 0.2
rho ~= 1 to 2 g cm^-3
```

That is a warning region, not a hard boundary.  A generated table can contain
points there, but it should not be accepted for production without checking
for branch dependence and thermodynamic pathologies.

## Practical extension policy

1. Rebuild a pilot grid with upstream FreeEOS 3.0.0 over the current MESA
   FreeEOS table ranges first.

2. Confirm that the MESA table-format quantities agree with the existing
   2.2.1-derived tables to the expected small level away from ionization
   features.  The upstream release notes say 3.0.0 and 2.2.1 should give
   essentially the same results when 2.2.1 converged; most differences should
   come from constants, ionization potentials, and bug fixes.

3. Extend into OPAL-like coverage and overlap deliberately.  The existing
   MESA OPAL/SCVH controls use:

   ```text
   logT_all_OPAL = 7.5
   logT_all_HELM = 7.6
   logQ_min_OPAL_SCVH = -10.09
   logQ_max_OPAL_SCVH = 5.3
   logRho1_OPAL_SCVH_limit = 3.50
   ```

   FreeEOS is a plausible table source in this overlap, especially for
   non-OPAL compositions, but it should still be compared against OPAL/SCVH,
   HELM, PC, or Skye in the regions where those EOSes are expected to be more
   reliable.

4. Treat deeper high-density, low-temperature pressure ionization as
   research-only until validated.  The first target could be a limited
   low-mass-star / hot-brown-dwarf extension, not a blanket replacement for
   HELM/PC/Skye.

5. Reject or mask any region showing:

   ```text
   increasing-density and decreasing-density scans land on different branches
   increasing-temperature and decreasing-temperature scans land on different branches
   discontinuous species fractions or EOS derivatives
   negative or inconsistent heat capacities
   non-smooth Gamma1, chiRho, chiT, grad_ad
   loss of Maxwell consistency from the table interpolant
   ```

## Low-temperature low-density FreeEOS edge

The point

```text
T = 100 K
log10(T/[K]) = 2.0
log10(rho/[g cm^-3]) = -15
```

has:

```text
log10(Q) = log10(rho) - 2*log10(T) + 12
         = -15 - 2*2 + 12
         = -7
```

so it is inside the nominal `logQ` span of the existing FreeEOS table layout.
There is no density/pressure-ionization argument against covering this point
with FreeEOS.  It is extremely dilute, so the non-ideal terms should be tiny
and the EOS should approach an ideal neutral or molecular gas limit.

The real FreeEOS questions are:

1. Does upstream FreeEOS 3.0.0 converge cleanly with the chosen option suite
   at `log10(T)=2.0`, `log10(rho)=-15`, and the desired composition?
2. Are the returned derivatives smooth and thermodynamically consistent?
3. Does the low-temperature chemistry included by FreeEOS cover the physics
   needed for the use case?

The third point matters because FreeEOS is primarily a stellar-interior EOS.
It includes the 20 supported elements and important H species such as `H2` and
`H2+`, but it is not a general cool-atmosphere chemistry package with an
extended molecular/dust network.  At `T = 100 K`, a solar or metal-rich gas
may need molecules not represented by FreeEOS.  For pure H/He or for a simple
ideal/molecular low-density tail, FreeEOS may be adequate if it converges.

The local MESA tree does not appear to contain a separate cool-chemistry EOS
table set that would supersede FreeEOS in this corner.  The cool-chemistry
tables found locally are opacity tables in `kap`, for example Ferguson,
Freedman, C/O-enhanced, and AESOPUS machinery:

```text
kap/preprocessor/src/ferguson.f90
kap/preprocessor/src/freedman.f90
kap/private/kap_aesopus.f90
kap/AESOPUS_AGSS09.h5
```

Those are opacity data, not thermodynamic EOS data.  On the EOS side, the
local low-temperature molecular information is mainly from H/He EOS sources
such as SCVH and from FreeEOS itself.  For example, `scvh_core.f90` reads H2
number concentration from the SCVH-derived tables, but this is not a broad
cool molecular chemistry network:

```text
eos/eosDT_builder/src/scvh_core.f90
```

Therefore, within this repository, the choices for `T ~= 100 K` and very low
density are essentially:

```text
FreeEOS if it converges and its included species are adequate
SCVH/H-He style molecular information for H/He only where available
an explicit ideal neutral/molecular EOS tail
a new external cool-chemistry EOS source, if broader molecules are required
```

There is one historical caution: an old FreeEOS mailing-list report from
2008 found an EOS1 low-temperature convergence issue near `log10(T)=3.0115`
for hydrogen-bearing mixtures over a broad density range, and noted that the
`H2` and `H2+` partition-function data were known to be good down to
`log10(T)=3`.  FreeEOS 3.0.0 has many reliability fixes, so this old report is
not a hard limit for current FreeEOS, but it is enough reason to require a
direct 3.0.0 point/probe scan before promising clean tables at `log10(T)=2`.

Thus, for a pure FreeEOS table project, the correct statement is:

```text
T = 100 K, log10(rho) = -15 is not excluded by the rho,Q geometry or by
pressure-ionization physics.  It is a low-temperature convergence and
chemistry-validity test for FreeEOS 3.0.0.
```

If FreeEOS 3.0.0 converges and passes derivative checks there, the table can
cover it.  Whether the physics is "better" than a deliberately constructed
ideal/molecular low-temperature tail depends on the required chemistry.

## Other public low-temperature EOS resources

There are public EOS or chemistry resources near this general subject, but no
obvious ready-made table found so far has the exact combination wanted here:

```text
T ~= 100 K
rho ~= 1d-15 g cm^-3
stellar-evolution thermodynamic quantities and derivatives
arbitrary H/He/metal compositions
MESA-style rho,T table integration
```

Relevant public resources found:

- CMS19 / Chabrier, Mazevet, & Soubiran H/He EOS tables are public and cover
  `T = 1d2..1d8 K`, but their density range is quoted as
  `rho = 1d-8..1d6 g cm^-3`.  That misses `rho = 1d-15 g cm^-3` by seven
  dex and is H/He, not a general metal mixture.
- REOS.3 H and He tables are public as paper supplemental material and cover
  `T = 60..1d7 K`, but the quoted density range is
  `rho = 1d-10..1d3 g cm^-3`.  That still misses `rho = 1d-15 g cm^-3` by
  five dex and is H/He.
- AESOPUS is relevant to low-temperature gas, arbitrary mixtures, chemistry,
  and opacities, but its commonly quoted range is `3.2 <= log10(T) <= 4.5`.
  It does not cover `T = 100 K`, and the MESA-local AESOPUS usage is opacity,
  not an EOS replacement table.
- GGchem and FastChem are public equilibrium-chemistry solvers that reach
  `T = 100 K` and include much broader molecule/condensate chemistry than
  FreeEOS.  They are not drop-in thermodynamic EOS tables; they are potential
  ingredients for constructing an ideal-gas/molecular EOS tail.
- Warm-dense-matter or HED EOS databases such as FPEOS, SESAME-like tables,
  and core-collapse nuclear EOS tables are not targeted at this low-density,
  cold stellar-mixture corner.  Their temperature and density ranges or
  material scope do not match the requested region.

Thus the current working interpretation is: public tables exist for H/He
planetary/interior EOS work, and public chemistry solvers exist for cool gas,
but a public, broad-composition, low-density, `T ~= 100 K`, MESA-ready
thermodynamic EOS table is not obvious.  That makes FreeEOS 3.0.0 plus direct
validation, or a constructed ideal/molecular tail using chemistry-solver
outputs, the plausible paths.

## Actual chemistry path

If the goal is actual cool chemistry rather than a pre-existing EOS table,
then FreeEOS is not the only possible source.  Public equilibrium chemistry
codes such as GGchem and FastChem are more chemically complete in this
temperature regime than FreeEOS.  For example, GGchem is described as working
down to `100 K` with many elements, molecules, and condensates; later
atmosphere-grid work describes GGchem as applicable over `100..6000 K`.
FastChem / PyFastChem reports tests over:

```text
100 K <= T <= 6000 K
1d-13 bar <= P <= 1d3 bar
```

For the specific low-density point:

```text
T = 100 K
rho = 1d-15 g cm^-3
P = rho k_B T / (mu m_u)
```

For a molecular-gas value `mu ~= 2.3`,

```text
P ~= 3.6d-6 dyn cm^-2 ~= 3.6d-12 bar
```

so the point is above the quoted FastChem pressure floor.  The pressure
depends on `mu`, which is part of the chemical-equilibrium solution, but this
order-of-magnitude check says the target point is not obviously outside the
chemistry-solver pressure domain.

The catch is that these are chemistry solvers, not complete MESA EOS tables.
They can provide equilibrium species abundances and mean molecular weight
under an ideal-gas assumption.  To make an EOS suitable for MESA-style use,
one would still need to construct thermodynamic quantities and derivatives.

At very low density, such as `rho = 1d-15 g cm^-3`, this ideal chemical EOS is
physically well motivated because non-ideal interactions are negligible:

```text
P = n_tot k_B T = rho k_B T / (mu m_u)
```

where `mu = mu(T,rho,{X_i})` is determined by chemical equilibrium.  The
energy and entropy cannot be just a constant-gamma ideal gas if molecules,
ionization, dissociation, and condensates are active.  They need the species
partition functions / thermochemical potentials:

```text
f(T,rho,{X_i}) -> P, S, E, Cv, chiRho, chiT, Gamma1
```

The clean version would again be a free-energy construction: use the chemical
equilibrium solution to evaluate a mixture Helmholtz or Gibbs free energy,
interpolate that as the primitive, and derive all EOS quantities from it.  A
less clean but practical pilot path would tabulate chemistry outputs over
`rho,T`, finite-difference them carefully, and verify Maxwell/derivative
consistency.

Therefore, if actual cool chemistry is required, the likely path is not to
search for a hidden complete EOS table.  It is to build a low-density
ideal-chemical EOS table using GGchem/FastChem-style equilibrium chemistry
and a free-energy-based table formulation.

The working generation pipeline should be:

```text
chemistry solver
  -> species fractions and mu(T,rho,{X_i})
thermochemistry / free energy
  -> P, E, S and derivatives from one potential
table generator
  -> MESA-compatible EOS table
```

The transition plan would be:

```text
cool chemistry EOS:       ~100 K to ~6000 K, low-density ideal/molecular gas
FreeEOS / OPAL-like EOS:  hotter stellar-envelope/interior regime
HELM/PC/Skye-like EOS:    fully ionized, high-T/high-density plasma regime
```

The exact handoff temperature should be validated from overlapping chemistry,
ionization fraction, and derivative smoothness rather than hard-coded from
the nominal solver range.

## FastChem docs and inputs

Primary FastChem documentation:

```text
https://newstrangeworlds.github.io/FastChem/
```

The useful documentation sections for a MESA EOS table builder are:

```text
Installation
  -> Obtaining the code
  -> Prerequisites
  -> Configuration and compilation with CMake
  -> Installation of pyFastChem with Python

Standard input and output files
  -> Element abundance file
  -> Species data files
  -> Basic element data file
  -> FastChem parameter file
  -> Output files

Running pyFastChem
  -> Provided Python examples
  -> Detailed steps for running pyFastChem
  -> Output functions of pyFastChem

Detailed pyFastChem module description
  -> pyFastChem constructor
  -> pyFastChem input and output structures
  -> pyFastChem functions

Internal parameters
  -> Standard parameters
  -> Advanced parameters
```

Availability notes:

```text
FastChem docs:   https://newstrangeworlds.github.io/FastChem/
FastChem GitHub: https://github.com/NewStrangeWorlds/FastChem
PyPI package:    https://pypi.org/project/pyfastchem/
```

As of the current check, PyPI lists `pyfastchem 3.1.3`, released
2025-07-31.  The docs describe FastChem Cond / FastChem 3.x and include both
C++ standalone and Python workflows.  For this table project, the Python
workflow is the right first driver because it lets the table builder manage
grid traversal, pressure iteration, caches, diagnostics, and validation
without modifying MESA or compiling MESA.

FastChem inputs that must be controlled for reproducibility:

```text
element abundance file
gas species thermochemical data files
condensate species data files, if using FastChem Cond
basic element data file
FastChem parameter file
condensation/rainout mode
solver tolerances and fallback settings
```

The table builder should copy or checksum all of these into the generated
table metadata so the EOS table is reproducible.

## Chemistry solver choice

Working recommendation: start the table-builder prototype with FastChem Cond
/ PyFastChem, and use GGchem as a validation/completeness cross-check.

Reasoning:

- FastChem Cond is the easiest starting point for a table-building workflow.
  It has a current Python package (`pyfastchem 3.1.3`, released 2025-07-31),
  C++ source, CMake, documented Python examples, and a current GitHub
  repository.  The Python layer is useful for quickly generating `rho,T`
  grids, composition sweeps, caches, and diagnostics before committing to a
  Fortran/C++ integration.
- FastChem is designed for speed.  FastChem 2 reports being up to 50 times
  faster than FastChem 1 over the tested grid, and FastChem 3.1 adds a
  multidimensional Newton fallback for harder gas-phase cases.  This matters
  for large table generation.
- FastChem Cond now includes condensation and rainout options, with about
  290 condensate species added in the Cond paper.  This is enough for a first
  cool molecular/condensation EOS-tail experiment.
- GGchem appears chemically richer.  The ASCL record describes up to 40
  elements, up to 1155 molecules, and up to 200 condensates from NIST-JANAF
  and SUPCRTBL.  The 2018 paper emphasizes robust low-temperature
  thermochemical data down to `100 K`, and later atmosphere work describes
  GGchem as a Gibbs-free-energy-based equilibrium chemistry tool over
  `100..6000 K`.
- GGchem is Fortran 90, which is attractive for eventual MESA-adjacent
  integration, but its packaging/API surface is less immediately convenient
  than PyFastChem for a first table generator.

Choice by goal:

```text
fastest first prototype / easiest availability:
  FastChem Cond via PyFastChem

production table generator after the method is proven:
  FastChem C++ library or standalone driver, possibly wrapped by a small
  reproducible table-build script

richest chemistry/data cross-check:
  GGchem

Fortran-native integration experiment:
  GGchem

thermodynamic consistency:
  neither solver alone is sufficient; consistency comes from constructing
  one free-energy / thermochemistry primitive and deriving P, E, S, and
  derivatives from it
```

The last point is important.  FastChem and GGchem solve chemical equilibrium.
They do not automatically provide a complete MESA EOS with `P, E, S, Cv,
chiRho, chiT, Gamma1` and consistent derivatives.  The table project must own
that thermodynamic layer.  Solver choice mainly determines species fractions,
mean molecular weight, condensates, available thermochemical data, speed, and
ease of driving large grids.

The practical plan should be:

1. Use PyFastChem to generate a small pilot `rho,T` grid for the intended
   MESA mixtures.
2. Build the ideal-chemical thermodynamic layer around those species
   abundances and thermochemical data.
3. Run derivative/Maxwell checks and compare against FreeEOS/SCVH in overlap
   regions.
4. Repeat the same grid with GGchem where possible.  If GGchem and FastChem
   disagree materially in regions important to MESA, inspect species lists,
   condensates, element depletion assumptions, and thermochemical data before
   choosing one for production.

### What GGchem has that FastChem does not

GGchem's advantages for this project are mostly chemistry breadth and
`rho,T` convenience:

- GGchem is advertised as handling up to 40 elements, up to 1155 molecules,
  and up to 200 condensates from NIST-JANAF and SUPCRTBL.  This is the
  richer published species set compared with FastChem's standard setup.
- GGchem can be run with either pressure or mass density constrained.  That is
  directly useful for MESA-style `rho,T` table generation.  FastChem is
  naturally a `P,T` chemistry solver, so a `rho,T` table would need an outer
  solve for pressure because `mu(P,T,{X_i})` feeds back into
  `rho = P mu m_u / (k_B T)`.
- GGchem has several selectable equilibrium-constant data files and can
  combine them, with later files overriding earlier ones.  That makes it
  useful for sensitivity tests of the thermochemical data.
- GGchem outputs detailed condensate diagnostics such as supersaturation,
  condensed-unit concentrations, remaining gas-phase element abundances,
  dust-to-gas mass ratio, and dust volume per H nucleus.
- GGchem uses a Fortran code base, which may be convenient for an eventual
  MESA-adjacent integration after a table workflow is proven.

FastChem's advantages are mostly engineering and robustness for table
production:

- FastChem Cond has an actively packaged C++/Python implementation,
  `pyfastchem` on PyPI, CMake support, documentation, examples, and a current
  release series.  The GitHub page lists FastChem 3.1.3 as the latest release
  on 2025-07-31.
- FastChem 3.x includes equilibrium condensation and rainout, automatic
  stable-condensate selection, and about 290 liquid/solid condensates in the
  Cond paper.
- FastChem is specifically designed around speed.  It has a semi-analytic gas
  solver, an updated multidimensional Newton fallback in 3.1, and a Python
  interface that makes grid generation and failure diagnostics easier.
- FastChem is tested over `100..6000 K` and `1d-13..1d3 bar` for solar
  abundances.

Do we need GGchem?

```text
For the first table-builder prototype: no, probably not.
For validation and chemistry uncertainty: yes, very useful.
For final production if FastChem misses species important to the target
mixtures or fails in rho,T corners: possibly.
```

The minimum sensible path is:

1. Prototype with PyFastChem because it is easier to install, drive, and
   automate.
2. Build the `rho,T` outer pressure solve around FastChem and check whether
   it is stable enough.
3. Run GGchem on selected grid cuts, especially `T < 1000 K`, high C/O,
   metal-rich, and condensation-heavy mixtures.
4. Promote GGchem from cross-check to production only if it materially
   changes `mu`, major species, condensate depletion, or derivative behavior
   in regions that MESA would actually use.

## Condensates versus rock EOS

The chemistry solvers can include condensates, but this is not the same as
switching to a condensed-matter or rock equation of state.

In FastChem Cond or GGchem, condensates are chemical-equilibrium species.
They affect:

```text
gas-phase element depletion
mean molecular weight of the remaining gas
latent heat / thermochemical energy if included in the EOS construction
opacity/cloud inputs if coupled to opacity or atmosphere calculations
```

They do not by themselves provide:

```text
solid/liquid pressure-density relation for rock/ice/metal
elastic/compressional properties of condensed phases
high-pressure mineral phase transitions
planetary-interior solid/liquid EOS behavior
```

For the very low-density stellar-atmosphere corner, condensates should mostly
be treated as chemically removing material from the gas and contributing
thermochemical terms if the table includes them.  The gas pressure is still
the ideal pressure of the gas-phase species:

```text
P_gas = n_gas k_B T
```

If a zone becomes dominated by condensed material and the condensed phase
needs to carry pressure or internal energy like a solid/liquid, then a
separate material EOS is required.  That is a different problem from adding a
cool chemical gas tail.  Candidate external sources would be planetary /
materials EOS tables for rock, ice, iron, or mixtures, not FastChem/GGchem
alone.

## MESA-format table build control flow

The easiest first target is not a new MESA runtime EOS.  It is a generator
that emits a text table in the same shape as current `eosDT` / FreeEOS tables,
then lets the existing MESA loader and bicubic interpolation machinery read
it.  This is the lowest-friction path for testing.

Current MESA `eosDT` / FreeEOS text-table shape:

```text
header line with column names
version X Z num_logTs logT_min logT_max del_logT num_logQs logQ_min logQ_max del_logQ

for each logQ:
  blank/subheader lines
  logQ value
  row header
  rows containing:
    logT
    logPgas
    logE
    logS
    chiRho
    chiT
    Cp
    Cv
    dE_dRho
    dS_dT
    dS_dRho
    mu
    log10_free_e
    gamma1
    gamma3
    grad_ad
    eta
    optional extra diagnostic columns ignored by current loader
```

The existing reader only loads the first 16 values after `logT`.  Current
FreeEOS tables append diagnostics such as `MESA`, `logRho`, `dpe`, `dsp`, and
`dse`; the loader ignores them.  A FastChem table can do the same for
diagnostics, but any value needed by MESA at runtime must be in the first 16
loaded quantities or require a code change.

Relevant local code:

```text
eos/private/eosdt_load_tables.f90
  jlogPgas ... jeta define the 16 loaded file columns
  Load1_eosDT_Table reads the text table
  Make_XEoS_Interpolation_Data converts log10 file quantities to MESA internals

eos/eosFreeEOS_builder/src/free_eos_table.f90
  write_table shows the current FreeEOS table writer pattern
```

Important unit convention:

```text
file logPgas, logE, logS, log10_free_e are base-10 logs
MESA internal values for lnPgas, lnE, lnS, lnfree_e are natural logs
the loader multiplies those table logs by ln(10)
```

Control flow for an easy prototype:

```text
1. Choose composition grid.
   Use MESA-like X,Z points first, then later add richer mixture parameters
   if needed.

2. Choose rho,T grid.
   MESA table files are indexed by logQ and logT, where
   logQ = logRho - 2*logT + 12.
   For each logQ,logT point recover:
   logRho = logQ + 2*logT - 12.

3. For each grid point, solve chemistry.
   FastChem is naturally P,T, so for a rho,T table:

     guess P
     run FastChem(P,T,composition)
     get gas species, condensates, mu
     rho_calc = P * mu * m_u / (k_B*T)
     iterate P until rho_calc matches target rho

   Use previous neighboring grid points as pressure guesses to keep the
   sweep fast and continuous.

4. Compute the thermodynamic primitive.
   Preferred path:

     evaluate a mixture free energy f(T,rho,{X_i})
     include translational, internal partition-function, dissociation,
     ionization, and condensate thermochemical terms as appropriate

   The chemistry solver supplies equilibrium composition; the EOS layer owns
   the free-energy / thermodynamic construction.

5. Derive MESA quantities from one primitive.
   Use automatic differentiation or high-order controlled finite differences
   on the free-energy interpolant to obtain:

     P, S, E, Cv, chiRho, chiT, gamma1, gamma3, grad_ad

   Then write the 16 MESA file columns in the expected order.

6. Append diagnostics.
   Extra columns can include:

     free_energy
     logRho
     pressure-solve iterations
     FastChem status code
     max element-conservation residual
     dominant gas species
     dominant condensates
     ionization/electron fraction summaries
     Maxwell residuals dpe,dsp,dse

   These are for validation unless the MESA loader is extended.

7. Load with existing MESA table machinery.
   For the first prototype, mimic the FreeEOS builder output and install into
   a separate data directory or a controlled local branch.  Do not overwrite
   shipped data until the validation suite is convincing.

8. Validate before runtime use.
   Check smoothness in logQ/logT, Maxwell identities, positive heat
   capacities, gamma behavior, species continuity, element conservation, and
   overlap agreement with SCVH/FreeEOS/HELM where their physics should apply.
```

For performance, the table build should be split into deterministic stages:

```text
chemistry cache:
  raw FastChem results at each grid point

thermo cache:
  free-energy primitive and derived MESA quantities

text table:
  MESA-compatible table plus diagnostics

validation report:
  residuals, failure masks, species maps, handoff diagnostics
```

This makes it possible to rerun the expensive chemistry only when species
inputs or solver settings change.

## Ionization and species tabulation

Ionization states and species fractions are worth saving, but they should not
be confused with the thermodynamic path.

For a free-energy EOS, the runtime thermodynamic outputs should come from one
thermodynamic primitive:

```text
f(T,rho,{X_i}) -> P, E, S, Cv, chiRho, chiT, gamma1, gamma3, grad_ad
```

The species fractions, ionization fractions, electron density, and condensate
fractions are equilibrium coordinates that explain the minimum.  They are
useful for diagnostics and possibly for coupling to opacities or rates, but
they should not be independently interpolated to derive thermodynamics unless
the derivative consistency is handled carefully.

Recommended output policy:

```text
main MESA table:
  only the 16 quantities currently loaded by eosdt_load_tables.f90

auxiliary diagnostic table:
  species number fractions
  element partitioning between gas and condensates
  electron fraction / free electron density
  dominant ionization states
  dominant molecules and condensates
  chemistry solver status and residuals

future MESA extension:
  add a formal auxiliary species table only if runtime code needs these
  quantities, e.g. for opacity coupling or atmosphere diagnostics
```

For `log10_free_e`, the existing MESA EOS table format already includes a
free-electron quantity.  A FastChem-based low-temperature table should compute
that consistently from the chemistry solution.  In mostly neutral molecular
gas this value can become extremely small, so it may need the same kind of
flooring/diagnostic treatment used elsewhere in MESA to avoid interpolation
pathologies.

## Coupling to extended FreeEOS

The intended long-term architecture should probably be a two-source EOS table
family:

```text
FastChem chemical EOS:
  cool, low-density, molecule/condensate-rich gas

extended FreeEOS:
  warm/hot stellar-envelope and OPAL/SCVH-like domain

HELM/PC/Skye:
  high-temperature, fully ionized, dense plasma regimes
```

The FastChem table should not try to replace FreeEOS in the whole OPAL/SCVH
region.  Conversely, extended FreeEOS should not be forced to solve all the
way into the cold `T ~= 100 K` chemistry corner if a dedicated chemical EOS is
available.  The clean plan is to give them an overlap and blend only where
both are demonstrably valid.

The natural overlap target is roughly:

```text
FastChem tested range:  100 K <= T <= 6000 K
FreeEOS useful low-T edge after validation: around log10(T) ~= 3.0..3.8+
candidate overlap:     log10(T) ~= 3.2..3.7, low-density gas
```

The actual overlap should be chosen from validation, not from nominal solver
limits.  Useful handoff diagnostics:

```text
species state:
  H2/H/H+/e- fractions
  major molecules such as CO, H2O, TiO, etc.
  condensate depletion fractions

thermodynamics:
  P, E, S
  chiRho, chiT
  Cv, Cp
  gamma1, gamma3, grad_ad
  Maxwell residuals

numerics:
  chemistry convergence status
  FreeEOS convergence status
  derivative smoothness in logT and logRho/logQ
```

If both FastChem and FreeEOS tables are built from free-energy primitives, the
best handoff is a free-energy-level blend in an overlap band:

```text
f_blend = w(T,rho) * f_FastChem + (1 - w(T,rho)) * f_FreeEOS
```

where `w` is a smooth function with zero first derivative at the ends of the
blend region.  All thermodynamic quantities should then be derived from
`f_blend`, including the derivative terms from the blend function itself.
This avoids blending `P`, `E`, `S`, and `Gamma1` independently.

If only MESA-style derived-column tables are available at first, a temporary
derived-quantity blend can be used for testing, but it should be treated as a
prototype path because independent blending of derived quantities does not
guarantee Maxwell consistency.

Implementation sequence:

```text
1. Build FastChem free-energy table over the cool chemistry range.
2. Build regenerated FreeEOS table over the current FreeEOS/OPAL-like range.
3. Validate both independently.
4. Define an overlap mask where both are valid.
5. Compare f, P, E, S, and derivatives in the overlap.
6. Choose smooth blend weights in logT/logRho or logT/logQ.
7. Derive final MESA columns from the blended free-energy table.
8. Emit MESA-compatible table plus auxiliary species diagnostics.
```

This keeps the FastChem work connected to the FreeEOS-extension work without
making either solver responsible for physics outside its strongest regime.

## Validation equations

If a free-energy table is used, the interpolation should satisfy at least the
following checks:

```text
P = rho^2 * (df/drho)_T
S = - (df/dT)_rho
E = f + T*S
Cv = (dE/dT)_rho = T * (dS/dT)_rho
```

and a Maxwell relation:

```text
(dS/drho)_T = - (d/drho)(df/dT)_rho
```

The pressure derivative consistency check can be written:

```text
chiRho = (d ln P / d ln rho)_T
chiT   = (d ln P / d ln T)_rho
```

The existing table generator already computes the MESA-format quantities
from FreeEOS outputs in `free_eos_eval`; a free-energy-based generator would
instead use these equations as the primary path and use FreeEOS returned
derivatives as validation or fallback.

## Code references

- `docs/source/eos/overview.rst`: documents MESA FreeEOS as based on
  FreeEOS 2.2.1, EOS1, and a domain similar to OPAL.
- `eos/eosFreeEOS_builder/README`: describes the current MESA FreeEOS table
  builder and notes that it falls back to MESA EOS if FreeEOS fails.
- `eos/eosFreeEOS_builder/inlist_solar`: current builder grid,
  `log10T = 2.1..8.2` and `log10Q = -10.09..5.69`.
- `eos/eosFreeEOS_builder/src/free_eos_table.f90`: current MESA table
  generator.  `free_eos_eval` calls upstream `free_eos`, and `write_table`
  writes the tabulated MESA quantities.
- `eos/private/eosdt_load_tables.f90`: current text-table reader.  The
  `jlogPgas` through `jeta` constants define the 16 values loaded from each
  row, and `Make_XEoS_Interpolation_Data` converts base-10 table logs to
  MESA natural-log internals.
- `eos/defaults/eos.defaults`: runtime FreeEOS selection and blending
  controls, including the current FreeEOS table upper `logT` and cut/blend
  controls.
- `eos/private/eosdt_eval.f90`: `Get_FreeEOS_alfa` applies the runtime
  FreeEOS region and cut/blend logic.
- `eos/public/eos_def.f90`: current MESA FreeEOS `X,Z` grid metadata.

## Source references

- MESA EOS documentation:
  `https://docs.mesastar.org/en/stable/eos/overview.html`
- FreeEOS file releases:
  `https://sourceforge.net/projects/freeeos/files/freeeos/`
- FreeEOS 3.0.0 release notes:
  `https://sourceforge.net/projects/freeeos/files/freeeos/3.0.0%20Source/`
- Irwin Paper V on convergence and pressure-ionization issues:
  `https://freeeos.sourceforge.net/convergence.pdf`
- CMS19 H/He EOS abstract and public-table statement:
  `https://ore.exeter.ac.uk/repository/handle/10871/36246`
- REOS.3 H/He EOS abstract and supplemental-table statement:
  `https://www.osti.gov/biblio/22340126`
- AESOPUS low-temperature gas opacity/EOS paper:
  `https://www.aanda.org/articles/aa/full_html/2009/48/aa12598-09/aa12598-09.html`
- GGchem code record:
  `https://www.ascl.net/2104.018`
- GGchem A&A paper:
  `https://www.aanda.org/articles/aa/full_html/2018/06/aa32193-17/aa32193-17.html`
- MSG atmosphere-grid paper describing GGchem `100..6000 K` applicability:
  `https://www.aanda.org/articles/aa/full_html/2024/10/aa50108-24/aa50108-24.html`
- FastChem paper:
  `https://academic.oup.com/mnras/article/479/1/865/5035838`
- FastChem 2 paper with `100..6000 K`, `1d-13..1d3 bar` testing:
  `https://academic.oup.com/mnras/article/517/3/4070/6702739`
- FastChem Cond paper:
  `https://academic.oup.com/mnras/article/527/3/7263/7424171`
- FastChem documentation:
  `https://newstrangeworlds.github.io/FastChem/`
- FastChem GitHub repository:
  `https://github.com/NewStrangeWorlds/FastChem`
- PyFastChem package notes with the same broad test range:
  `https://pypi.org/project/pyfastchem/`

## Open tasks

- [ ] Inspect PyFastChem examples and API details for the exact output fields
      needed by the chemistry cache.
- [ ] Decide the first FastChem species/condensate data set and checksum all
      input data files for reproducibility.
- [ ] Implement a standalone `rho,T -> P,T -> FastChem -> mu` pressure solve
      without modifying MESA.
- [ ] Decide the free-energy / thermochemistry primitive and derivative
      method for the first prototype.
- [ ] Emit a tiny MESA-shaped text table with the 16 required EOS columns and
      appended diagnostics.
- [ ] Load the table through existing MESA EOS table code only after explicit
      permission to compile/run.
- [ ] Compare the FastChem table against SCVH/FreeEOS in overlap regions.
- [ ] Run selected GGchem cross-checks before promoting FastChem chemistry to
      production.
