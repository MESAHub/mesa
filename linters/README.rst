=======
Linters
=======


MESA Linters
============

These scripts help check MESA for common problems.
The linters must be run from ``$MESA_DIR/linters``.

check_photos.py
---------------

This script makes sure that all of the values included in
``star_data_step_input.inc`` are saved in photos.

check_columns.py
----------------

This script checks that the history and profile code and the
``*_columns.list`` files are in sync. It can also check the column files
in the test_suite if enabled in the source code.

check_defaults.py
-----------------

This script checks that the controls and star_job defaults, the
``*_controls.inc``, and ``*_io.f90`` files are in sync.

check_pgstar.py
----------------

This script checks that the inlist_pgstar files in the test_suite
have valid history/profile columns. Note it will have false positives for things only
defined in a run_star_extras.f90 and does not check non inlist_pgstar files
pgstar sections.


fix_underlines.py
-----------------

This script checks that the ~-level underlines in the defaults files
are the correct lengths and fixes any problems it finds.

mesa_linter.py
--------------

Checks fortran files for consistency with MESA's style guide.

Run over files and provide details of failures

.. code-block:: console

  python mesa-linter.py *.f90

Provides only a count per file of the number of issues found

.. code-block:: console

  python mesa-linter.py -s *.f90

check_test_suite_onwers.py
--------------------------

This script checks that each test suite case is listed in CODEOWNERS and
checks whether each test case has 0 owners (very bad) or just 1 owner
(less bad but should be fixed).

list_test_owner.py
------------------

Parses CODEOWNERS and prints the test cases owned by each person.

check_stop.py
-------------

Checks all .f90 files for stop 1 or stop 'string' and replaces them with a call
to mesa_error(__FILE__,__LINE__) . This way we can always find where an error
occurs rather than the non-unique stop 1 location.

check_empty_writes.py
---------------------

Checks all .f90 files for writes that are both empty and unformatted (i.e write (*,*) ) (it also handles if the write is to a unit).
Unformatted writes are non-portable and can cause issues with ifort.

check_omp_critical.py
---------------------

Checks all .f90 files for unnamed critical blocks. When critical blocks are unnamed they acts as one block
thus each must be run separately. When critical blocks are named they can each be run in parallel improving performance.

update_columns.py
-----------------

Copies the history and profile default columns file into each test case while preserving enabled
options in each test case

fix_inlists.py
--------------

Fixes various controls in the test_suite inlists that should not be enabled by default.

update_ctrls.py
---------------

This handles the replacement of ``s%`` with ``s% ctrl%`` while being smart about only doing that to controls.


Code Linters
============

We use a variety of linters to check the MESA codebase for style and
correctness. These linters are run automatically via GitHub actions
on pushes and pull requests, and will fail if there are any issues.

Fortran linter
--------------

We use `fortitude <https://github.com/PlasmaFAIR/fortitude>`__ to
lint the Fortran files in the MESA codebase for style and correctness.
The linter can be run in the commandline with

.. code-block:: console

  fortitude check

Python linter
-------------

We use `ruff <https://docs.astral.sh/ruff/>`__ to lint and check
Python files.
The linter can be run in the commandline with

.. code-block:: console

  ruff format --check
  ruff check

reStructuredText linter
-----------------------

We use `sphinx-lint <https://pypi.org/project/sphinx-lint/>`__ to lint
reStructuredText files.
The linter can be run in the commandline with

.. code-block:: console

  sphinx-lint --ignore=gyre/gyre
