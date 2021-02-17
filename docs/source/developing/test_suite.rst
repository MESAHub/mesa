**************
The Test Suite
**************

The MESA test suites live in ``star/test_suite``, 
``binary/test_suite``, and ``astero/test_suite``.

Running Tests
-------------

You can count the total number of available tests::

  ./count_tests

You can list the available tests::

  ./list_tests

You can get the name of a test (or tests) by specifying one or more integers as arguments::

  ./list_tests <N>

You can run all tests::

  ./each_test_run

You can run an individual test by specifying a single integer, corresponding to the number from the list of tests::

  ./each_test_run <N>

After a test runs, the file ``out.txt`` will contain the concatenated
output from stdout and stderr.


Anatomy of a test
-----------------

In the test case directory 
^^^^^^^^^^^^^^^^^^^^^^^^^^

All MESA test suite cases must be designed to be multi-part---that is
to run multiple inlists---even if they only run a single inlist in
practice.

The script that runs a single part is called ``rn1`` (which usually
looks like the ``rn`` script in the stock work directory).  The
restart script ``re`` is usually stock and so restarts a single part.

The ``rn`` script is responsible for the logic related to running
multiple parts.  Generally, each test part has a header inlist
(``inlist_*_header``), which then points to other inlists. At the
start of each part, this will be copied to ``inlist``, and will be the
file that MESA reads as its primary inlist.


A good test should be able to regenerate its starting model.  When
then environment variable ``MESA_RUN_OPTIONAL`` is set, test cases
should regenerate their initial model.  The ``rn`` script should look
something like

.. literalinclude:: ../../../star/test_suite/test_case_template/rn
   :language: bash

It is essential to make use of ``test_suite_helpers`` as this ensures
that important information is produced in a TestHub-friendly format.

A test case must not write to stderr.  The presence of output on
stderr will cause the test to be classified as a failure.  This
restriction catches ``stop`` statements, calls to ``mesa_error``, or
other error conditions.

A good test should run relatively quickly.  Costly parts can be
skipped over using a saved model and only run when
``MESA_RUN_OPTIONAL`` is set.

A good test should check its stopping condition before producing an
output model (i.e., set ``required_termination_code_string``).

A good test should use ``run_star_extras`` to perform more detailed
physical and/or numerical checks or report longitudinally interesting
values to the TestHub (see :ref:`MESA Test Hub`).

Some tests should only be run in certain circumstances (e.g., if GYRE
is installed, if OP MONO opacities are present).  Such a test should
still compile and run even when these conditions are not met, but
should output the string "this test was intentionally skipped".  When
this string is present in the test output, the test scripts will
consider the test a success and no further checks will be performed.


In the test suite directory
^^^^^^^^^^^^^^^^^^^^^^^^^^^

The existence of a sub-directory in ``test_suite`` is not sufficient
to tell MESA to perform a given test.  The master list of tests lives
in the file ``test_suite/do1_test_source``.  This has an a series of
entries like::

  do_one 1.3M_ms_high_Z "stop because log_surface_luminosity >= log_L_upper_limit" "final.mod" x300

The 4 arguments to ``do_one`` are:

#. test name (this should be the directory name)
#. termination condition (this string must exist in the test output for the test to be considered a success)
#. model filename (this is a file from the end of the run that will have its checksum computed)
#. photo filename (MESA will restart from this file and check that the checksum of the file was the same after the run and the restart)

The model filename argument permits the special value ``skip`` which causes this check to be skipped.
The photo filename argument permits the special value ``skip`` which causes this check to be skipped, or ``auto`` which restarts from the antepenultimate file (determined by filesystem modification times).

Once an entry in ``do1_test_source`` has been added, the test can be
run as described in `Running Tests`_.

Documenting a test
------------------

Files documenting the test suite exist in two places.  Look to
existing test cases for detailed examples.  (Currently
:ref:`conductive_flame` is set up this way.)

In the test case directory
^^^^^^^^^^^^^^^^^^^^^^^^^^

.. include:: ../../../star/test_suite/test_case_template/README.rst


In the docs directory
^^^^^^^^^^^^^^^^^^^^^

Once the ``README.rst`` file is created, a file in
``docs/source/test_suite`` should be created.  This is the location of
the include directive that incorporates the contents of the
``README.rst`` file.  This also includes the anchor to the test suite
for use in cross-references.

An entry linking to the test case page should be included in
``docs/source/test_suite.rst``.  This page will eventually list all
cases.  Briefly describing the test cases all in one place will make
it easier for people to find a test of interest.
    
The description should be one sentence broadly describing what the
case does and then one optional sentence about anything interesting
illustrated by the case (e.g., other_hooks).
   
MESA Test Hub
-------------

The `MESA Test Hub <https://testhub.mesastar.org/>`__ records the
results of nightly test runs.  The owner/maintainer of the MESA Test
Hub is Bill Wolf.

When a test is run through ``each_test_run``, a file ``testhub.yml``
will be produced.  This is the information that will be reported to
the TestHub by ``mesa_test``.  The file will look similar to this:

.. code-block:: none

    ---
    test_case: make_co_wd
    module: :star
    omp_num_threads: 2
    inlists:
        - inlist: inlist_co_core_header
          runtime_minutes:   2.29
          steps:    151
          retries:      0
          redos:      0
          solver_calls_made:    151
          solver_calls_failed:      0
          solver_iterations:    846
          log_rel_run_E_err:        -3.8240036219155993
        - inlist: inlist_remove_env_header
          runtime_minutes:   2.08
          steps:    438
          retries:     18
          redos:      0
          solver_calls_made:    456
          solver_calls_failed:     18
          solver_iterations:   1361
          log_rel_run_E_err:        -3.8220395900315363
        - inlist: inlist_settle_header
          runtime_minutes:   3.55
          steps:    161
          retries:      0
          redos:      0
          solver_calls_made:    161
          solver_calls_failed:      0
          solver_iterations:    735
          log_rel_run_E_err:        -7.2371589938451679
    mem_rn: 7814204
    success_type: :run_test_string
    mem_re: 3471908
    success_type: :photo_checksum
    checksum: 48f31b30ef89bb579de45f68919d57be
    outcome: :pass

The output is collected in a variety of places.  The highest level
information (i.e. no indent) that summarizes the test case itself
comes from ``each_test_run``.  The inlist name information (i.e. lines
with ``-``) comes from ``do_one`` used in the test suite ``rn``
scripts (which is provided by ``test_suite_helpers``).

.. note ::

    If a test suite case contains a file named ``.ignore_checksum``,
    the checksum value will not be reported to the TestHub.

The per-inlist information about the performance of MESA is provided
by the ``test_suite_extras`` (see
``star/job/test_suite_extras_def.inc`` and
``star/job/test_suite_extras.inc``) included in the
``run_star_extras`` / ``run_binary_extras`` of the test suite case.
In particular, in the star/binary ``after_evolve`` hook, the call to
``test_suite_after_evolve`` (and its callees) append the relevant info
to the ``testhub.yml`` file.

The aforementioned routines control output that should be present for
every MESA test case.  Some test cases may want to output additional
information.  To do so, set elements in the provided arrays
``testhub_extras_names`` and ``testhub_extras_values``.  The values in
these arrays (at the time the after evolve hook is called) will be
automatically added to ``testhub.yml`` and reported to the TestHub.

.. warning::

   The values of ``extra_testhub_names`` must be unique within each
   test case.  In a multi-part test, you cannot report a value with
   the same name in each part.  Use some extra identifier to break
   this ambiguity (i.e., ``he_core_mass_part1``,
   ``he_core_mass_part2``).


For example, in ``c13_pocket``, the ``run_star_extras`` sets::

    ! put target info in TestHub output
    testhub_extras_names(1) = 'max_c13'; testhub_extras_vals(1) = max_c13
    testhub_extras_names(2) = 'mass_max_c13' ; testhub_extras_vals(2) = mass_max_c13
    testhub_extras_names(3) = 'pocket_mass_c13'; testhub_extras_vals(3) = pocket_mass_c13
    testhub_extras_names(4) = 'delta_surface_c12'; testhub_extras_vals(4) = delta_surface_c12

which results in the additional output

.. code-block:: none

    - inlist: inlist_c13_pocket_header
      ...
      extra_testhub_names:
            - 'max_c13'
            - 'mass_max_c13'
            - 'pocket_mass_c13'
            - 'delta_surface_c12'
      extra_testhub_vals:
            -     7.5183958374071852E-02
            -     5.7873466406466934E-01
            -     3.3303926652486070E-05
            -     7.2015481753889745E-05

