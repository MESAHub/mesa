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

``./rn`` copies ``input/history.data`` to ``LOGS/history.data`` (unless that
file already exists), then runs the program. The program writes
``LOGS/post_proc_output.data`` and finishes with::

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
   the inlist sets it to ``.false.``, as a run made without inline colors
   usually does. To reuse the settings of the run that wrote the history file, point at
   that run's inlist instead of copying them::

      read_extra_colors_inlist(1) = .true.
      extra_colors_inlist_name(1) = '/path/to/that/run/inlist_colors'

``&post_proc``
   ``history_file``
      The history file to read. Default ``'LOGS/history.data'``.
   ``output_file``
      The file to write. Default ``'LOGS/post_proc_output.data'``.

The history file must contain these columns:

* ``star_age``
* ``log_Teff`` or ``Teff``
* ``log_g``
* ``log_R`` or ``radius``
* ``surface_h1`` and ``surface_he4``

The default ``history_columns.list`` does not include ``surface h1`` and
``surface he4``. Add them to the run's history columns before running it.

The surface composition sets the atmosphere metallicity. The program maps
:math:`Z/X` onto the atmosphere-table metallicity axis as
:math:`\log_{10}((Z/X)/z\_over\_x\_ref)`, with :math:`Z = 1 - X - Y`, and
clamps the result to the range of the table. This is the same mapping, and
the same ``z_over_x_ref`` control, that ``star`` uses for the inline colors
history columns.

Output
======

``output_file`` has one header line of column names, then one row per history
row. The first four columns are ``star_age``, ``Teff``, ``log_g`` and
``m_div_h`` (the metallicity passed to the atmosphere interpolation). Then
come the colors history columns: ``Mag_bol``, ``Flux_bol``, ``Interp_rad``
and one magnitude per filter.

Post-processing your own run
============================

#. Add ``surface h1`` and ``surface he4`` to the run's history columns and
   run it.
#. Either set ``history_file`` in ``&post_proc`` to that run's
   ``LOGS/history.data``, or copy it in with::

      ./grab_post_inputs.sh /path/to/that/run

   which copies the run's ``inlist*`` files here and its history file to
   ``LOGS/history.data``, after backing up any files it replaces. The run's
   ``inlist`` then supplies ``&colors``, and since it has no ``&post_proc``
   the defaults apply.
#. ``./rn``

Plotting
========

``python_helpers/stitch_history.py`` appends the colors columns to the
history file named in ``&post_proc``. It renames the original file with an
``_old`` suffix. After stitching, the plotting scripts in
``star/test_suite/custom_colors/python_helpers`` work on the result. Run them
from ``colors_post/python_helpers``, since they read ``../LOGS/history.data``.

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
