************
Common Tasks
************

This page describes commonly-encountered development tasks.

Inlist Options
--------------

Adding a new option
^^^^^^^^^^^^^^^^^^^

.. note::

   The simplest approach is often to find a similar existing option
   and follow the pattern revealed by searching for the occurrences of
   that option in the source code.

This guide applies to star-related options (``star_job``, ``controls``, and ``pgstar``).

First, add an entry of the relevant type in the appropriate ``.inc``
file in ``star_data/private``.

* ``star_job`` : ``star_job_controls.inc``
* ``controls`` : ``star_controls.inc``
* ``pgstar`` : ``pgstar_controls.inc``

Next, add an entry with the default option value and some
documentation (see :ref:`reference/format:Format for MESA defaults
files`) into the appropriate ``*.defaults`` file in ``star/defaults``.

* ``star_job`` : ``star_job.defaults``
* ``controls`` : ``controls.defaults``
* ``pgstar`` : ``pgstar.defaults``

Finally, update the appropriate ``*_io.f90`` files in
``star/private``, which are responsible for reading/writing the Fortran
namelists.

* ``star_job`` : ``star_job_ctrls_io.f90``
* ``controls`` : ``ctrls_io.f90``
* ``pgstar`` : ``pgstar_ctrls_io.f90``

This will require several edits.

* Add your option to the ``namelist`` declaration near the start of that file.

* Add code transferring your option from the value read by the namelist to the ``star_info`` structure in the ``store_*`` routine.

* For the case of ``star_job`` and ``controls``, add code transferring
  your option from the ``star_info`` structure to the namelist value
  in the ``set_*_for_writing`` routine.  This is what is used to
  enable the options :ref:`reference/star_job:save_star_job_namelist`
  and :ref:`reference/controls:save_controls_namelist`.


History/Profile Output
----------------------

.. note::

   The simplest approach is often to find a similar existing output
   and follow the pattern revealed by searching for the occurrences of
   that quantity in the source tree.

Adding a history column
^^^^^^^^^^^^^^^^^^^^^^^

First, update ``star/private/star_history_def.f90``.  Here, you need
to define the integer parameter of the form ``h_my_new_column``
internally used to identify the column.  Later in the file, assign the
name it will be given via an assignment like
``history_column_name(h_my_new_column) = 'my_new_column'``.

Next, update ``star/private/history.f90`` to add a new case in the
switch statement, ``case(h_my_new_column)``.  The code in this case
should set the variable ``val`` to the desired output value.g

Finally, add the new column to ``star/defaults/history_columns.list``
along with a short comment about what this column means.  It should
likely be commented out (off by default).


Adding a profile column
^^^^^^^^^^^^^^^^^^^^^^^

First, update ``star/private/star_profile_def.f90``.  Here, you need
to define the integer parameter of the form ``p_my_new_column``
internally used to identify the column.  Later in the file, assign the
name it will be given via an assignment like
``profile_column_name(p_my_new_column) = 'my_new_column'``.

Next, update ``star/private/profile_getval.f90`` to add a new case in
the switch statement, ``case(p_my_new_column)``.  The code in this
case should set the variable ``val`` to the desired output value
corresponding to cell ``k``.

Finally, add the new column to ``star/defaults/profile_columns.list``
along with a short comment about what this column means.  It should
likely be commented out (off by default).
