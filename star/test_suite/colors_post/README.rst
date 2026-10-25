.. _colors_post:

***********
colors_post
***********

.. tags:: star, colors

This test suite case computes synthetic photometry for an existing MESA
history file, without rerunning the stellar evolution. It is useful when a run
was made with ``use_colors = .false.``, or to recompute magnitudes for a
different filter set, atmosphere grid, distance or magnitude system.

Unlike the other test suite cases, ``colors_post`` does not build ``star``. It
builds a small standalone program, ``src/run.f90``, that reads a history file
and calls the ``colors`` module directly for every row.

Running the test
================

.. code-block:: bash

   make clean
   make
   ./rn

``./rn`` reads ``input/history.data`` and writes ``LOGS/history.data``, a
copy of it with the colors columns added. It finishes with::

   colors_post: successfully post-processed all N models

``input/history.data`` is a trimmed history from the ``custom_colors`` test
case (a 7 |Msun| star from the pre-main sequence to core helium depletion),
run with ``surface h1`` and ``surface he4`` added to its history columns.

Inputs
======

Both namelists are read from ``inlist``.

``&colors``
   The same controls as for a run with inline colors; see
   ``colors/defaults/colors.defaults``. ``use_colors`` is turned on even if
   the inlist leaves it ``.false.``, as a run made without inline colors
   usually does.

``&post_proc``
   ``history_file``
      The history file to read. Required.
   ``output_file``
      The file to write. Default ``'LOGS/history.data'``.
   ``fallback_m_div_h``
      The [M/H] to use for every row when the history file has no
      ``surface_h1`` and ``surface_he4`` columns. Ignored when it has them.
      Not set by default, in which case such a history file is an error.

The history file must contain ``log_Teff`` (or ``Teff``), ``log_g`` and
``log_R`` (or ``radius``). It should also contain ``surface_h1`` and
``surface_he4``, which the default ``history_columns.list`` does not include;
add ``surface h1`` and ``surface he4`` to the run's history columns, or set
``fallback_m_div_h``. If ``model_number`` is present it is passed on to the
colors module, which uses it to name per-model SED files.

The surface composition sets the atmosphere metallicity. The program maps
:math:`Z/X` onto the atmosphere-table metallicity axis as
:math:`\log_{10}((Z/X)/z\_over\_x\_ref)`, with :math:`Z = 1 - X - Y`, and
clamps the result to the range of the table. This is the same mapping, and
the same ``z_over_x_ref`` control, that ``star`` uses for the inline colors
history columns. ``fallback_m_div_h`` is clamped the same way.

Output
======

``output_file`` is a MESA history file: the header and every column of
``history_file``, plus the colors history columns ``Mag_bol``, ``Flux_bol``,
``Interp_rad`` and one magnitude per filter. A colors column that already
exists in ``history_file``, for example from a run with inline colors, is
replaced rather than duplicated, so running again on the output gives the
same file.

The output is written to a temporary file that replaces ``output_file`` only
once every row is done. A failed run leaves any earlier output untouched, and
``output_file`` may be the same file as ``history_file``.

Post-processing your own run
============================

#. Copy this directory and run ``make`` there.
#. In ``inlist``, set ``history_file`` to the run's ``LOGS/history.data``.
#. Set the ``&colors`` controls, or reuse the run's::

      read_extra_colors_inlist(1) = .true.
      extra_colors_inlist_name(1) = '/path/to/that/run/inlist_colors'

   Extra inlist names are resolved relative to the current directory, so
   name the file that contains the run's ``&colors`` namelist itself, not a
   top-level ``inlist`` that refers to it by a relative name.
#. ``./rn``

Plotting
========

The plotting scripts in ``star/test_suite/custom_colors/python_helpers`` read
``../LOGS/history.data``, so they work on the output when run from
``colors_post/python_helpers``.

``python_helpers/plot_zero_points.py`` runs the post-processing once each for
the Vega, AB and ST magnitude systems and overlays the resulting CMDs. It
changes ``mag_system`` in ``./inlist`` only, so ``mag_system`` must be set
there rather than in an extra inlist.

Differences from inline colors
==============================

Gravity and radius are the same in both: the history ``log_g`` and ``log_R``
columns are the outermost-cell values that ``star`` passes to the colors
module. Two inputs differ, so post-processed magnitudes can differ from
inline ones:

* ``star`` passes the temperature of the outermost cell, ``T(1)``.
  ``colors_post`` passes ``Teff``.
* ``star`` takes :math:`X` and :math:`Z` from the photosphere cell.
  ``colors_post`` uses the surface (outermost cell) abundances.
