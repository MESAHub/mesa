colors unit test
================

Compiles and runs a standalone Fortran program (``src/test_colors.f90``) that
exercises the MESA colors module without requiring a full stellar evolution run.
The test uses the Kurucz2003 atmosphere grid, the Johnson filter set, and the
Vega photometric system.


Running
-------

::

    ./mk          # compile
    ./rn          # run, writing output to stdout
    ./ck          # run and diff against test_output

``ck`` delegates to ``utils/test/ck``, which pipes stdout to ``tmp.txt`` and
compares it against ``test_output`` with ``diff -b``.

Regenerating test_output
------------------------

Run the test on Linux and redirect stdout::

    ./rn > test_output

The reference file must always be generated on Linux.  See the note below on
why.

Output format and platform precision
--------------------------------------

All numerical output uses the Fortran format ``1pe26.16`` (16 significant
figures, full double precision).  This is intentional: it allows future
changes to the interpolation or integration routines to be detected even when
the perturbation is very small.

However, IEEE 754 arithmetic does not guarantee bit-identical results across
compiler/OS combinations when the computation involves accumulated
floating-point operations.  In practice, macOS and Linux can disagree in the
final one or two significant figures — differences of order 10\ :sup:`−15`
relative.  Because the shared ``ck`` script uses plain ``diff -b`` (exact
string comparison), any such disagreement causes the test to fail on macOS
even when the result is numerically correct.

The ``test_output`` reference file is therefore generated on Linux, and the
CI macOS job is expected to be more fragile at this precision level.  This is
a known limitation of using exact string comparison with 16-digit output.

Diagnostic plot
---------------

``python_helpers/plot_unit_test.py`` reads ``test_output`` and produces a
four-panel PGSTAR-styled PDF (``unit_test_diagnostic_pgstar.pdf``):

- colour indices (U−B, B−V, V−I, V−J) vs [M/H], log g, and Teff
- solar SED sanity panel with Johnson filter transmissions, a scaled
  Planck function, and the Wien-law peak

Run from ``colors/test/``::

    python3 python_helpers/plot_unit_test.py

Requires matplotlib and numpy.  If ``$MESA_DIR`` is set, the Johnson filter
curves are loaded from the colors data directory and overlaid on the SED panel;
otherwise the SED is plotted alone.