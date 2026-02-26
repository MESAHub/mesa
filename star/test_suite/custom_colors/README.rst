.. _custom_colors:

******
Colors
******

This test suite case demonstrates the functionality of the MESA ``colors`` module.

Running the Test Suite
======================


The ``custom_colors`` test suite case evolves a 7 M☉ star from the pre-main sequence through core hydrogen burning, up to the point of X_He,c < 0.01. 

The test is not scientifically rigorous—the pre-MS relaxation settings are aggressive and the mesh is fine—its purpose is solely to exercise the colors module across a wide range of stellar parameters (T_eff, log g) in a reasonable wall-clock time.


The ``custom_colors`` test suite is a standard MESA work directory. Before running it for the first time—or after making changes to the ``colors`` module source—the binary must be compiled.

``make clean``
   Removes all previously compiled object files and binaries from the ``build/`` directory. Run this before recompiling to ensure a clean state, particularly after switching MESA versions or modifying source files.

``make``
   Compiles the test suite and links it against the installed MESA libraries (including ``libcolors``). When it completes successfully it produces the ``build/bin/star`` executable.

``./rn``
   Runs the compiled stellar evolution model. MESA evolves the star according to the parameters in ``inlist_colors``, writing history and profile data to ``LOGS/`` and photometric outputs to ``SED/``.

A typical run workflow is:

.. code-block:: bash

   make clean   # wipe any previous build
   make         # compile
   ./rn         # run

If ``make`` completes without errors and ``./rn`` begins printing timestep output, the colors module is working correctly.



What is MESA colors?
====================

MESA colors is a runtime module that allows users to generate observer-ready data directly from stellar evolution models. Instead of limiting output to theoretical quantities like Luminosity (L) and Surface Temperature (T_eff), the colors module computes:

* **Bolometric Magnitude** (M_bol)
* **Bolometric Flux** (F_bol)
* **Synthetic Magnitudes** in specific photometric filters (e.g., Johnson V, Gaia G, 2MASS J).

How does the MESA colors module work?
=====================================

1.  At each timestep, the module takes the star's current surface parameters—Effective Temperature (T_eff), Surface Gravity (log g), and Metallicity ([M/H])—and queries a user-specified library of stellar atmospheres (defined in ``stellar_atm``). It interpolates within this grid to construct a specific Spectral Energy Distribution (SED) for the star's current parameters.

2.  This SED is then convolved with filter transmission curves (defined in ``instrument``) to calculate the flux passing through each filter.

3.  The fluxes are converted into magnitudes using the user-selected magnitude system (AB, ST, or Vega).

Inlist Options & Parameters
===========================

The colors module is controlled via the ``&colors`` namelist. Below is a detailed guide to the key parameters.

use_colors
----------

**Default:** ``.false.``

Master switch for the module. Must be set to ``.true.`` to enable any photometric output.

**Example:**

.. code-block:: fortran

   use_colors = .true.


instrument
----------

**Default:** ``'data/colors_data/filters/Generic/Johnson'``

Path to the filter instrument directory, structured as ``facility/instrument``.

* The directory must contain an index file with the same name as the instrument
  (e.g., ``Johnson``), listing one filter filename per line.
* The module loads every ``.dat`` transmission curve listed in that index and
  creates a corresponding history column for each.

.. rubric:: Note on paths

All path parameters (``instrument``, ``stellar_atm``, ``vega_sed``) are resolved
using the same logic:

* ``'data/colors_data/...'`` — no leading slash; ``$MESA_DIR`` is prepended. This
  is the recommended form for all standard data paths.
* ``'/absolute/path/...'`` — tested on disk first; if found, used as-is. If not
  found, ``$MESA_DIR`` is prepended (preserves backwards compatibility).
* ``'./local/path/...'`` or ``'../up/one/...'`` — used exactly as supplied,
  relative to the MESA working directory.

**Example:**

.. code-block:: fortran

   instrument = 'data/colors_data/filters/GAIA/GAIA'


stellar_atm
-----------

**Default:** ``'data/colors_data/stellar_models/Kurucz2003all/'``

Path to the directory containing the grid of stellar atmosphere models. Paths may be relative to ``$MESA_DIR``, relative to the working directory, or absolute. This directory must contain:

1.  **lookup_table.csv**: A map linking filenames to physical parameters (T_eff`, log g, [M/H]).
2.  **SED files**: The actual spectra (text or binary format).
3.  **flux_cube.bin**: (Optional but recommended) A binary cube for rapid interpolation.

The module queries this grid using the star's current parameters. If the star evolves outside the grid boundaries, the module will clamp to the nearest edge.

**Example:**

.. code-block:: fortran

   stellar_atm = 'data/colors_data/stellar_models/sg-SPHINX/'


distance
--------

**Default:** ``3.0857d19`` (10 parsecs in cm)

The distance to the star in centimetres, used to convert surface flux to observed flux.

* **Default Behaviour:** At 10 parsecs (3.0857 * 10^19 cm) the output is **Absolute Magnitudes**.
* **Custom Usage:** Set this to a specific source distance to calculate Apparent Magnitudes.

**Example:**

.. code-block:: fortran

   distance = 5.1839d20


make_csv
--------

**Default:** ``.false.``

If set to ``.true.``, the module exports the full calculated SED at every profile interval.

* **Destination:** Files are saved to the directory defined by ``colors_results_directory``.
* **Format:** CSV files containing Wavelength vs. Flux.
* **Use Case:** Useful for debugging or plotting the full spectrum of the star at a specific evolutionary age.

**Example:**

.. code-block:: fortran

   make_csv = .true.


colors_results_directory
------------------------

**Default:** ``'SED'``

The folder where CSV files (if ``make_csv = .true.``) and other outputs are saved.

**Example:**

.. code-block:: fortran

   colors_results_directory = 'sed'


mag_system
----------

**Default:** ``'Vega'``

Defines the zero-point system for magnitude calculations. Options are:

* ``'AB'``: Based on a flat spectral flux density of 3631 Jy.
* ``'ST'``: Based on a flat spectral flux density per unit wavelength.
* ``'Vega'``: Calibrated such that Vega has magnitude 0 in all bands.

**Example:**

.. code-block:: fortran

   mag_system = 'AB'


vega_sed
--------

**Default:** ``'data/colors_data/stellar_models/vega_flam.csv'``

Required only if ``mag_system = 'Vega'``. Points to the reference SED file for Vega, used to compute photometric zero-points. Paths may be relative to ``$MESA_DIR``, relative to the working directory, or absolute.

**Example:**

.. code-block:: fortran

   vega_sed = '/path/to/my/vega_SED.csv'


Data Preparation (SED_Tools)
============================

The ``colors`` module requires pre-processed stellar atmospheres and filter
profiles organised in a specific directory structure. To automate this
entire workflow, we provide the dedicated repository:

**Repository:** `SED_Tools <https://github.com/nialljmiller/SED_Tools>`_

SED_Tools downloads, validates, and converts raw spectral atmosphere grids and
filter transmission curves from the following public archives:

* `SVO Filter Profile Service <http://svo2.cab.inta-csic.es/theory/fps/>`_
* `MAST BOSZ Stellar Atmosphere Library <https://archive.stsci.edu/prepds/bosz/>`_
* `MSG / Townsend Atmosphere Grids <https://www.astro.wisc.edu/~townsend/msg/>`_

These sources provide heterogeneous formats and file organisations. SED_Tools
standardises them into the exact structure required by MESA:

* ``lookup_table.csv``
* Raw SED files (text and/or HDF5)
* ``flux_cube.bin`` (binary cube for fast interpolation)
* Filter index files and ``*.dat`` transmission curves

SED_Tools produces:

.. code-block:: text

    data/
    ├── stellar_models/<ModelName>/
    │   ├── flux_cube.bin
    │   ├── lookup_table.csv
    │   ├── *.txt / *.h5
    └── filters/<Facility>/<Instrument>/
        ├── *.dat
        └── <Instrument>

These directories can be copied or symlinked directly into your MESA
installation.

A browsable online mirror of the processed SED_Tools output is also available:

`SED Tools Web Interface (mirror) <https://nillmill.ddns.net/sed_tools/>`_

This server provides a live view of:

* All downloaded stellar atmosphere grids
* All available filter facilities and instruments
* File counts, disk usage, and metadata
* Direct links to the directory structure used by MESA


Defaults Reference
==================

Below are the default values for all user-facing ``colors`` module parameters as defined in ``colors.defaults``.

.. code-block:: fortran

      use_colors = .false.
      instrument = 'data/colors_data/filters/Generic/Johnson'
      stellar_atm = 'data/colors_data/stellar_models/Kurucz2003all/'
      vega_sed = 'data/colors_data/stellar_models/vega_flam.csv'
      distance = 3.0857d19  ! 10 parsecs in cm (Absolute Magnitude)
      make_csv = .false.
      sed_per_model = .false.
      colors_results_directory = 'SED'
      mag_system = 'Vega'

Visual Summary of Data Flow
===========================

.. code-block:: text

   +----------------+
   |   MESA Model   |
   | (Teff, logg, Z)|
   +----------------+
           |
           v
   +-------------------------------------------------------------------------+
   |                        MESA COLORS MODULE                               |
   | 1. Query Stellar Atmosphere Grid with input model                       |
   | 2. Interpolate grid to construct specific SED                           |
   | 3. Convolve SED with filters to generate band flux                      |
   | 4. Apply distance flux dilution to generate bolometric flux -> Flux_bol |
   | 5. Apply zero point (Vega/AB/ST) to generate magnitudes                 |
   |                                    (Both bolometric and per filter)     |
   +-------------------------------------------------------------------------+
           |
           v
   +----------------------+
   |    history.data      |
   | age, Teff, ...       |
   | Mag_bol, Flux_bol    |
   | V, B, I, ...         |
   +----------------------+


=====================
=====================
Python Helper Scripts
=====================

All Python helpers live in ``python_helpers/`` and are run from that directory.
All paths default to ``../LOGS/`` and ``../SED/`` relative to that location, matching the standard test suite layout.

plot_history_live.py
--------------------

**Purpose:** Live-updating four-panel diagnostic viewer for ``history.data``, designed to run *during* a MESA simulation.

**What it shows:**

* **Top-left**: Color–magnitude diagram (CMD) constructed automatically from whichever filters are present in the history file. Filter priority for color index selection follows: Gaia (Gbp−Grp), then Johnson (B−R or B−V), then Sloan (g−r), with a fallback to the first and last available filter.
* **Top-right**: Classical HR diagram (Teff vs. log L).
* **Bottom-left**: Color index as a function of stellar age.
* **Bottom-right**: Light curves for all available filter bands simultaneously.

All four panels are color-coded by MESA's ``phase_of_evolution`` integer if that column is present in the history file. If it is absent, points are colored by stellar age using a compressed inferno colormap that emphasises the most recent evolution.

**How to use:**

.. code-block:: bash

   cd python_helpers
   python plot_history_live.py

The script polls ``../LOGS/history.data`` every 0.1 seconds and updates the plot whenever the file changes. It will print a change notification for the first five updates, then go silent to avoid log spam. Close the window to exit.


plot_history.py
---------------

**Purpose:** Single-shot (non-live) version of the history viewer. Reads the completed ``history.data`` once and renders the same four-panel figure. Intended for post-run analysis and as a shared library imported by ``plot_history_live.py``, ``movie_history.py``, and ``movie_cmd_3d.py``.

**How to use:**

.. code-block:: bash

   cd python_helpers
   python plot_history.py

Produces a static figure and then calls ``plt.show()``. The script also exports ``MesaView``, ``read_header_columns``, and ``setup_hr_diagram_params`` for use by the other scripts.

**Note:** The first 5 model rows are skipped (``MesaView`` skip=5) to avoid noisy pre-MS relaxation artifacts at the very start of the run.

plot_sed_live.py
----------------

**Purpose:** Live-updating SED viewer that monitors the ``SED/`` directory for CSV files written by the colors module (requires ``make_csv = .true.`` in ``inlist_colors``).

**What it shows:** A single plot with a logarithmic wavelength axis showing:

* The full stellar SED (black line, plotted once from the first file found).
* The filter-convolved flux for each band (one colored line per filter).
* The Vega reference SED if a ``VEGA_*`` file is present.
* Colored background bands marking the X-ray, UV, optical, and IR regions of the electromagnetic spectrum.
* A text box in the corner showing the stellar mass, metallicity, and distance parsed from ``inlist_colors``.

The x-axis auto-scales to the wavelength range where the SED has flux above 1% of its peak; the y-axis scales to the range of convolved fluxes.

**How to use:**

.. code-block:: bash

   cd python_helpers
   python plot_sed_live.py

Runs live, refreshing every 0.1 s. Close the window or press Ctrl-C to stop. Can also be configured to save a video instead of displaying live by setting ``save_video=True`` in the ``SEDChecker`` constructor.


plot_sed.py
-----------

**Purpose:** Single-shot SED plot that reads all ``*SED.csv`` files in ``../SED/`` and overlays them in one figure. Simpler than ``plot_sed_live.py``—no live monitoring, no EM region shading, no inlist parsing.

**How to use:**

.. code-block:: bash

   cd python_helpers
   python plot_sed.py

Displays the combined SED figure. The x-axis is cropped to 0–60000 Å by default; edit the ``xlim`` argument in ``main()`` to change this.


plot_cmd_3d.py
--------------

**Purpose:** Interactive 3D scatter plot of the CMD with a user-selectable third axis. Useful for visualising how any history column correlates with the photometric evolution of the star.

**How to use:**

.. code-block:: bash

   cd python_helpers
   python plot_cmd_3d.py

On launch, the script prints all available columns from ``history.data`` and prompts you to type the name of the column to use for the Z axis (default: ``Interp_rad``, the interpolation radius in the atmosphere grid). Press Enter to accept the default or type any other column name. A rotatable 3D matplotlib window then opens showing color index vs. magnitude vs. your chosen column.


movie_history.py
----------------

**Purpose:** Renders the same four-panel display as ``plot_history_live.py`` into an MP4 video, with one frame per model row. Points accumulate from left to right in time, so the full evolutionary track builds up across the video.

**How to use:**

.. code-block:: bash

   cd python_helpers
   python movie_history.py

Writes ``history.mp4`` in the current directory at 24 fps, 150 dpi. Requires ``ffmpeg`` to be installed and accessible on ``$PATH``. Progress is shown with a ``tqdm`` progress bar if that package is available.


movie_cmd_3d.py
---------------

**Purpose:** Creates a 3D rotation video that starts looking straight down the ``Interp_rad`` axis (showing a plain CMD) and rotates to an oblique 3D perspective over 10 seconds, with 1-second holds at each end.

**How to use:**

.. code-block:: bash

   cd python_helpers
   python movie_cmd_3d.py

Writes ``cmd_interp_rad_rotation.mp4`` in the current directory at 30 fps. Requires ``ffmpeg``.


plot_newton_iter.py
-------------------

**Purpose:** Plots the per-Newton-iteration photometry data written to ``SED/iteration_colors.data`` by the colors module's solver-monitor hook. This file records the photometry at every Newton iteration within each timestep, capturing sub-timestep variability during convergence.

The script supports both an **interactive mode** (a terminal UI with a filterable column picker) and a **batch mode** driven by command-line arguments. Both modes support arbitrary column expressions (e.g., ``B-V``, ``(V-U)/(B-V)``, ``Teff/1000``) in addition to named columns. When a ``history.data`` file is present, the history track is overlaid on the plot in grey for comparison. Plots are saved as both PDF and JPG.

**Interactive mode:**

.. code-block:: bash

   cd python_helpers
   python plot_newton_iter.py

You will be prompted to select the plot type (2D scatter, 2D line, or 3D scatter), then the X, Y, (optionally Z), and color axes from a grid display. The picker supports substring filtering (``/text``), negative filtering (``!text``), and regex filtering (``//pattern``). You can also type a column name directly.

**Batch mode:**

.. code-block:: bash

   python plot_newton_iter.py -x iter -y V -c model
   python plot_newton_iter.py -x Teff -y "B-V" -c "U-B"
   python plot_newton_iter.py -x iter -y V -z R -c model       # 3D
   python plot_newton_iter.py --list-columns                   # show columns and exit

Key flags:

* ``-f`` / ``--file``: path to the iteration data file (default: ``SED/iteration_colors.data``)
* ``--history``: path to history file for overlay (default: ``../LOGS/history.data``)
* ``-x``, ``-y``, ``-z``, ``-c``: axis column names or expressions
* ``--cmap``: matplotlib colormap (default: ``viridis``)
* ``--no-show``: suppress the interactive window (save only)
* ``--list-columns``: print available columns and exit


movie_newton_iter.py
--------------------

**Purpose:** Creates an animated MP4 showing the Newton iteration data accumulating point by point over time. Iterations are sorted by model number then iteration number, so the animation follows the physical sequence of the simulation. History file points are overlaid incrementally—a new history point appears each time a model's full set of iterations is complete. Outliers are removed via iterative sigma clipping before rendering. Axis limits expand smoothly as new data arrives.

The script supports the same interactive and batch modes as ``plot_newton_iter.py``, and imports all shared functionality (terminal UI, data loading, expression parsing) directly from that module.

**Interactive mode:**

.. code-block:: bash

   cd python_helpers
   python movie_newton_iter.py

You will be prompted for X, Y, and color axes, video duration, FPS, output filename, and sigma-clipping threshold.

**Batch mode:**

.. code-block:: bash

   python movie_newton_iter.py -x Teff -y V -c model
   python movie_newton_iter.py -x Teff -y "B-V" -c iter --flip-y -d 60

Key flags:

* ``-f`` / ``--file``: path to iteration data file (default: ``SED/iteration_colors.data``)
* ``--history``: path to history file (default: ``../LOGS/history.data``)
* ``-x``, ``-y``, ``-c``: axis column names or expressions
* ``-d`` / ``--duration``: target video duration in seconds (default: 30)
* ``--fps``: frames per second (default: 30)
* ``--flip-y``: invert the Y axis (useful for magnitude axes)
* ``--sigma``: sigma-clipping threshold (default: 3.0)
* ``--cmap``: matplotlib colormap (default: ``viridis``)

Requires ``ffmpeg``. For GIF output, ``pillow`` can be used instead by changing the writer in the source.


plot_zero_points.py
-------------------

**Purpose:** Runs MESA three times in sequence with ``mag_system`` set to ``'Vega'``, ``'AB'``, and ``'ST'`` respectively, then overlays the resulting CMDs in a single comparison figure. Useful for quantifying the offsets between magnitude systems for a given set of filters.

**How to use:**

.. code-block:: bash

   cd python_helpers
   python plot_zero_points.py

The script must be run from within ``python_helpers/`` (it changes directory to ``../`` before calling ``./rn``). It temporarily modifies ``inlist_colors`` to set the magnitude system and to disable PGstar, restores the original file after each run, and saves each LOGS directory as ``LOGS_Vega``, ``LOGS_AB``, and ``LOGS_ST``. The comparison figure is saved to ``mag_system_comparison.png``.

.. warning::

   This script runs the full simulation three times. For large or slow models, consider reducing ``initial_mass`` or loosening ``varcontrol_target`` and ``mesh_delta_coeff`` inside the ``FAST_RUN_OVERRIDES`` dictionary at the top of the script before running.


run_batch.py
------------

**Purpose:** Runs MESA in batch over a list of parameter files, each of which defines a different stellar type (M dwarf, Pop III, OB star, etc.). After each run, the history file is renamed to avoid being overwritten by the next run.

**How to use:**

Edit ``param_options`` at the top of the script to list the ``extra_controls_inlist_name`` files you want to loop over, then run:

.. code-block:: bash

   cd python_helpers
   python run_batch.py

The script: (1) comments out all parameter entries in ``inlist_1.0``, (2) uncomments one entry at a time, (3) runs ``./clean``, ``./mk``, and ``./rn`` in the MESA work directory, (4) renames the resulting ``history.data`` to ``history_<param_name>.data``, and (5) re-comments all entries at the end. If a run fails, a warning is printed and the loop continues with the next parameter set.

**Note:** The MESA work directory is assumed to be one level above the location of ``run_batch.py`` (i.e., ``../``). The inlist to modify is ``inlist_1.0``. Adjust these paths at the top of the script if your layout differs.
