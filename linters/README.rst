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
*_columns.list files are in sync.


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


