=======
Linters
=======

These scripts help check MESA for common problems.

check_photos.py
---------------

This script makes sure that all of the values included in
``star_data_step_input.inc`` are saved in photos.


check_columns.py
----------------

This script checks that the history and profile code and the
``*_columns.list`` files are in sync. It can also check the column files
in the test_suite if enabled in the source code.

check_pgstar.py
----------------

This script checks that the inlist_pgstar files in the test_suite
have valid history/profile columns. Note it will have false postives for things only
defined in a run_star_extras.f90 and does not check non inlist_pgstar files
pgstar sections.


fix_underlines.py
-----------------

This script checks that the ~-level underlines in the defaults files
are the correct lengths and fixes any problems it finds.


mesa_linter.py
--------------

Checks fortran files for consistency with MESA's style guide.

Run over files and provide details of failures::

  python3 mesa-linter.py *.f90
	 
Provides only a count per file of the number of issues found::

  python3 mesa-linter.py -s *.f90


