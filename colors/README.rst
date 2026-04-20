.. _colors:

******
Colors
******

The MESA ``colors`` module calculates synthetic photometry and bolometric quantities during stellar evolution.

What is MESA colors?
====================

MESA colors is a runtime module that generates observer-frame photometry directly from a running stellar evolution model. At each timestep it interpolates a pre-computed stellar atmosphere grid to construct a Spectral Energy Distribution (SED) for the star's current parameters, convolves that SED with photometric filter transmission curves, and writes synthetic magnitudes and bolometric quantities to ``history.data``.

The outputs are:

* **Mag_bol** — bolometric magnitude, derived directly from the stellar luminosity
* **Flux_bol** — bolometric flux at the specified distance
* **Interp_rad** — distance in parameter space between the current stellar parameters and the nearest atmosphere grid point (diagnostic for interpolation quality)
* **One column per filter** — synthetic magnitude in every filter listed in the instrument index file, named by the filter filename (``*.dat`` suffix stripped)

This is fundamentally different from the pre-existing bolometric correction (BC) interface in MESA. The old BC approach interpolates a table of pre-computed magnitude offsets. The colors module instead constructs a full SED at the stellar parameters and performs the photometry in full—there is no intermediate bolometric correction step. The old BC interface functions (``get_bc_by_name``, ``get_abs_mag_by_id``, etc.) are retained as stubs in ``colors_lib.f90`` that return ``-99.9`` solely to satisfy the MESA interface; they are not called by the colors module itself.

How it works
============

Initialisation (once per run)
------------------------------

On startup, ``colors_setup_tables`` loads all data that will be needed at runtime:

1. **Atmosphere grid** — reads ``lookup_table.csv`` from the ``stellar_atm`` directory. This maps every SED filename in the grid to its (T_eff, log g, [M/H]) coordinates. Unique sorted grids and an O(1) ``grid_to_lu`` index map are built from this table.

2. **Flux cube** — attempts to load ``flux_cube.bin``, a pre-built 4D binary array of shape ``(n_T_eff, n_log_g, n_[M/H], n_λ)``. If the allocation succeeds, the entire atmosphere grid is held in RAM and all subsequent interpolation is done in memory. If the allocation fails (insufficient RAM), the module falls back to the per-file stencil path described below.

3. **Filter transmission curves** — loads every ``*.dat`` file listed in the instrument index file. Each filter's wavelength and transmission arrays are stored on the handle.

4. **Zero-points** — for each filter, all three zero-points (Vega, AB, ST) are computed once from the filter transmission and the Vega reference SED (``vega_sed``). The selected zero-point is then looked up at runtime rather than recomputed.

Per-timestep computation
-------------------------

At each history output step, ``data_for_colors_history_columns`` is called with the current stellar parameters: T_eff, log g, metallicity [M/H], radius R, and distance d.

**Step 1 — SED interpolation**

The module locates the containing cell in the (T_eff, log g, [M/H]) grid and interpolates to produce a flux array F_λ at the stellar surface. Two paths exist:

* **Flux cube path** (preferred): hermite tensor interpolation across the full pre-loaded 4D array. All lookups are in-memory array accesses.
* **Stencil fallback path** (low-RAM): an extended neighbourhood of SED files around the current grid cell is loaded on demand. Individual SED files are served from a bounded memory cache (256-slot circular buffer, ``sed_mem_cache_cap`` in ``colors_def.f90``) to avoid redundant disk reads. The stencil is invalidated and reloaded whenever the star moves into a new grid cell.

**Step 2 — Distance dilution**

The surface flux is scaled to the observer:

.. code-block:: text

   F_observed(λ) = F_surface(λ) × (R / d)²

where R is the stellar radius and d is ``distance`` (default 10 pc, giving absolute magnitudes).

**Step 3 — Bolometric quantities**

The bolometric flux is obtained by integrating the diluted SED over all wavelengths using adaptive Simpson's rule (falling back to the trapezoid rule for even-length arrays). The bolometric magnitude follows from the standard relation using the solar bolometric absolute magnitude.

**Step 4 — Synthetic photometry**

For each filter, the in-band flux is computed by integrating the product of the diluted SED and the filter transmission curve:

.. code-block:: text

   F_band = ∫ F_observed(λ) × T(λ) dλ  /  ∫ T(λ) dλ

The synthetic magnitude is then:

.. code-block:: text

   m = -2.5 × log10(F_band / F_zp)

where F_zp is the precomputed zero-point for the selected magnitude system (Vega, AB, or ST).

If the star's parameters fall outside the atmosphere grid, the module clamps to the nearest grid boundary. The ``Interp_rad`` column records the Euclidean distance (in normalised parameter space) between the stellar parameters and the nearest grid point.

Source files
============

.. code-block:: text

   colors/
   ├── public/
   │   ├── colors_def.f90       — data structures (Colors_General_Info, filter_data),
   │   │                          handle management, memory layout
   │   └── colors_lib.f90       — public API: colors_init, colors_shutdown,
   │                              alloc_colors_handle, colors_setup_tables,
   │                              history column interface, legacy BC stubs
   ├── private/
   │   ├── bolometric.f90       — bolometric magnitude and flux calculation
   │   ├── synthetic.f90        — per-filter convolution and magnitude calculation,
   │   │                          SED CSV output (make_csv / sed_per_model)
   │   ├── hermite_interp.f90   — hermite tensor interpolation (cube path)
   │   ├── linear_interp.f90    — trilinear interpolation (cube path fallback)
   │   ├── knn_interp.f90       — k-nearest-neighbour interpolation
   │   ├── colors_utils.f90     — I/O (SED, filter, lookup table, flux cube),
   │   │                          numerical integration, flux dilution,
   │   │                          stencil and SED memory cache management
   │   ├── colors_history.f90   — MESA history column interface
   │   └── colors_ctrls_io.f90  — &colors namelist I/O
   └── defaults/
       └── colors.defaults      — default values for all &colors parameters

File Formats
============

All file formats used by the module are plain text or simple unformatted binary, with no external library dependencies.

lookup_table.csv
----------------

A CSV file located at ``<stellar_atm>/lookup_table.csv``. One row per atmosphere model in the grid. Column names are read from the header and matched by name, so column order is flexible. The required columns are:

.. code-block:: text

   filename, Teff, logg, MH

where ``filename`` is the path to the SED file relative to the ``stellar_atm`` directory, ``Teff`` is in Kelvin, ``logg`` is log₁₀(g / cm s⁻²), and ``MH`` is [M/H].

SED files (atmosphere grid)
----------------------------

Plain two-column text files, one row per wavelength point:

.. code-block:: text

   wavelength(Å)   flux(erg/s/cm²/Å)

No header is required. All SED files within a given atmosphere grid must share the same wavelength grid; if they do not (e.g. BT-Settl), the stencil loader will interpolate non-conforming files onto the canonical wavelength grid of the first file loaded.

The colors data ships with the Kurucz2003 models for solar alpha and alpha = 0.4.

Spectral Grid Variants
~~~~~~~~~~~~~~~~~~~~~~

+-------------------------+--------+-----------+----------------------+-------------------+----------------+------------------------+
| Name                    | [α/Fe] | N Spectra | Teff (K)             | logg              | [M/H]          | Wavelength (Å)         |
+=========================+========+===========+======================+===================+================+========================+
| Kurucz2003all           | 0.0    | 3808      | 3500-50000 (76 pt)   | 0.00-5.00 (11 pt) | -2.50-0.50 (8 pt) | 147-1600000 (1199 pt) |
+-------------------------+--------+-----------+----------------------+-------------------+----------------+------------------------+
| Kurucz2003all__alpha_04 | 0.4    | 4284      | 3500-50000 (76 pt)   | 0.00-5.00 (11 pt) | -4.00-0.50 (9 pt) | 147-1600000 (1199 pt) |
+-------------------------+--------+-----------+----------------------+-------------------+----------------+------------------------+

**[α/Fe]** is the alpha-element enhancement - the abundance of O, Ne, Mg, Si, S, Ar, Ca, and Ti relative to iron compared to solar.
A value of 0.4 means those elements are enhanced by 0.4 dex above solar relative to iron.
The alpha-enhanced grid goes down to lower metallicity (-4.00 vs -2.50) because alpha-enriched stars are usually old, metal-poor Pop II stars.

`Castelli & Kurucz 2003 <https://arxiv.org/abs/astro-ph/0405087>`__

`SVO model <https://svo2.cab.inta-csic.es/theory/newov2/index.php?models=Kurucz2003all>`__

flux_cube.bin
-------------

An unformatted Fortran binary file produced by SED_Tools. It encodes the full atmosphere grid as a contiguous 4D array ``(n_Teff, n_logg, n_MH, n_λ)`` plus the corresponding axis coordinate arrays. Loading this file at startup eliminates all per-timestep disk I/O for the flux cube path. The file is optional; if absent, the module uses the stencil fallback path.

Filter files (.dat)
-------------------

Two-column text files, one row per wavelength point. Lines beginning with ``#`` are treated as comments and skipped:

.. code-block:: text

   # optional comment lines
   wavelength(Å)   transmission(0–1)

The instrument directory must also contain an index file whose name matches the instrument (e.g. ``Johnson``), listing one filter filename per line. The module loads every filter listed in this index.

Vega reference SED (vega_flam.csv)
------------------------------------

A two-column CSV with a single header line:

.. code-block:: text

   wavelength,flux
   ...

Wavelengths in Å, flux in erg/s/cm²/Å. Only the first two columns are read; additional columns are ignored. This file is required when ``mag_system = 'Vega'`` and is used solely to precompute filter zero-points at initialisation.

SED output files (make_csv)
----------------------------

When ``make_csv = .true.``, the module writes one CSV file per profile interval to ``colors_results_directory``. Each file contains the diluted (observer-frame) SED at that timestep:

.. code-block:: text

   wavelength(Å),flux(erg/s/cm²/Å)

If ``sed_per_model = .true.``, filenames are suffixed with the model number (e.g. ``SED_00042.csv``); otherwise the same filename is overwritten at each interval.

Data Preparation (SED_Tools)
============================

The ``colors`` module requires pre-processed stellar atmosphere grids and filter
profiles organised in a specific directory structure. The dedicated repository
`SED_Tools <https://github.com/nialljmiller/SED_Tools>`_ automates this entire
workflow.

SED_Tools downloads, validates, and converts raw spectral atmosphere grids and
filter transmission curves from the following public archives:

* `SVO Filter Profile Service <http://svo2.cab.inta-csic.es/theory/fps/>`_
* `MAST BOSZ Stellar Atmosphere Library <https://archive.stsci.edu/prepds/bosz/>`_
* `MSG / Townsend Atmosphere Grids <https://www.astro.wisc.edu/~townsend/msg/>`_

SED_Tools standardises these heterogeneous sources into the directory structure
and file formats required by MESA:

.. code-block:: text

    data/
    ├── stellar_models/<ModelName>/
    │   ├── lookup_table.csv
    │   ├── flux_cube.bin
    │   └── *.txt / *.h5
    └── filters/<Facility>/<Instrument>/
        ├── <Instrument>        ← index file
        └── *.dat

The resulting ``data/`` tree can be copied or symlinked directly into
``$MESA_DIR``, after which the default path values in ``colors.defaults`` will
resolve correctly without any further configuration.

A browsable mirror of the processed SED_Tools output is available at
`https://nillmill.ddns.net/sed_tools/ <https://nillmill.ddns.net/sed_tools/>`_,
providing a live view of all available atmosphere grids and filter facilities
along with file counts, disk usage, and metadata.
