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
*_columns.list files are in sync.  (NOTE: It currently exhibits some
false positives, especially related to RSP due to case-sensitivity.)

