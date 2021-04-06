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
    :start-line: 6


In the docs directory
^^^^^^^^^^^^^^^^^^^^^

Once the ``README.rst`` file is created, a link in ``docs/source/test_suite`` should be created.

.. code-block:: sh

  ln -s ../../../star/test_suite/test_case/README.rst test_case.rst

This instructs Sphinx to incorporate the contents of the ``README.rst`` file.

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
    omp_num_threads: 36
    run_optional: false
    fpe_checks: false
    inlists:
        - inlist: inlist_co_core_header
          runtime_minutes:   3.20
          model_number:         1008
          star_age:     4.3063267775397134E+08
          num_retries:            7
          log_rel_run_E_err:        -7.0353816429993685
          steps:         1008
          retries:            7
          redos:            0
          solver_calls_made:         1015
          solver_calls_failed:            7
          solver_iterations:        10280
        - inlist: inlist_remove_env_header
          runtime_minutes:   2.62
          model_number:         1568
          star_age:     4.3064418442484504E+08
          num_retries:           30
          log_rel_run_E_err:        -6.9640200618973518
          steps:          961
          retries:           34
          redos:            0
          solver_calls_made:          995
          solver_calls_failed:           34
          solver_iterations:         6403
        - inlist: inlist_settle_header
          runtime_minutes:   2.44
          model_number:         1746
          star_age:     4.3359174056506151E+08
          num_retries:           10
          log_rel_run_E_err:       -10.1236273392339484
          steps:          184
          retries:            4
          redos:            0
          solver_calls_made:          184
          solver_calls_failed:            0
          solver_iterations:         1097
    mem_rn: 11815904
    success_type: :run_test_string
    restart_photo: x650
    mem_re: 4635248
    success_type: :photo_checksum
    checksum: cb6df95a221722e7317a6e53c9c61272
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

This output is generated each time MESA terminates (except after a restart).  Therefore, the per-inlist quantities that can be reported by TestHub are those accessible within a single part of a MESA run.  By default, we report

   + ``runtime_minutes``
   + ``steps``
   + ``retries``
   + ``redos``
   + ``solver_calls_made``
   + ``solver_calls_failed``
   + ``solver_iterations``

In a multi-part test case, the per-part values can be summed to give the properties of the complete run.

TestHub also reports quantities that can reflect information preserved by MESA across parts.  These are transmitted via their inclusion in the model file.  That means the values reported by cases that use saved models to skip optional parts will be influenced by the performance at the time the saved model was generated.  Additionally, some parts may include inlist options that reset or modify these quantities.  TestHub reports the values at the end of each part, but the precise meaning of these quantities cannot be understood without reference to the details of the test case.

   + ``model_number``
   + ``star_age``
   + ``num_retries``
   + ``log_rel_run_E_err``

.. note ::

  These values can be useful when diagnosing test case issues because they directly correspond to quantities in the terminal output.

Some test cases may want to output additional
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

Continuous integration testing
------------------------------

Multiple developers have set up their machines to enable continuous integration testing. These machines will
automatically pull the changes in the repository, run the test suite, and report back to testhub.
To make more efficient the usage of these machines they will respond to certain keywords if found in the
commit message.

.. note ::
    It is up to each person providing the computing resources to implement each keyword. Thus some 
    machines will ignore these keywords and run the test suite normally. Therefore, these are only
    "requests" for the computing machines not "orders".

The message (with brackets) may appear anywhere in the commit message. 


[ci skip]
^^^^^^^^^

Compile ``MESA`` but do not run the test suite. Useful when changes only touch documentation or the changes can not affect the final result.

[ci split]
^^^^^^^^^^

Splits the running of the test suite between machines. Current, if set, cannon will run the first half of the test cases while helios will run the second half.

[ci optional]
^^^^^^^^^^^^^

Runs ``MESA`` with the environment variable ``MESA_RUN_OPTIONAL=t`` set. This requests that the slower optional parts of each test case be ran.

[ci optional n]
^^^^^^^^^^^^^^^

Where ``n`` is an integer. Same as ``[ci optional]`` but only run the first ``n`` test cases.
