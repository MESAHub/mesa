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

After a test runs, the files ``out.txt`` and ``err.txt`` will contain
the output from stdout and stderr respectively.


Philosophy of the tests
-----------------------

The MESA test suite serves multiple roles.  For developers, it
provides day-to-day checks that code changes did not cause regressions
and provides longer term opportunities to monitor the evolving
performance of the code.  For users, the test suite primarily serves
as a source of examples.

The coverage and quality of the test suite must be sufficient to
ensure that as MESA is developed, it retains its key capabilities
(e.g., the ability to evolve through the He flash or to robustly
evolve massive stars).  To this end, a number of test suite cases
descend from examples presented in MESA instrument papers.


Anatomy of a test
-----------------

.. note::

   If you want to add a new test case, the first step is to understand
   the layout and motivation of the standard tests, as described in
   this section.  As a starting point, a simple template exists at
   ``star/test_suite/test_case_template`` (and the example files
   rendered below are contained therein).  For more complex
   situations, you may want to take inspiration from existing test
   cases.


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


A good test should be able to regenerate its starting model.
The ``rn`` script should look something like

.. literalinclude:: ../../../star/test_suite/test_case_template/rn
   :language: bash

When the environment variable ``MESA_SKIP_OPTIONAL`` is set, some
parts of the test run may be skipped by copying standard versions
of saved models (which must be included in the test case).

It is essential to make use of ``test_suite_helpers`` as this ensures
that important information is produced in a TestHub-friendly format.

Similarly, it is essential that the ``run_star_extras`` in the test
suite case use the ``test_suite_extras`` (see
``star/job/test_suite_extras_def.inc`` and
``star/job/test_suite_extras.inc``).  The calls to these hooks (see below)
generate the necessary information for the :ref:`developing/test_suite:MESA Test Hub`.

.. literalinclude:: ../../../star/test_suite/test_case_template/src/run_star_extras.f90
   :language: fortran
   :start-at: subroutine extras_startup
   :end-at: end subroutine extras_startup
   :emphasize-lines: 9

.. literalinclude:: ../../../star/test_suite/test_case_template/src/run_star_extras.f90
   :language: fortran
   :start-at: subroutine extras_after_evolve
   :end-at: end subroutine extras_after_evolve
   :emphasize-lines: 9


Properties of a good test
~~~~~~~~~~~~~~~~~~~~~~~~~

A test case must not write to stderr.  The presence of output on
stderr will cause the test to be classified as a failure.  This
restriction catches ``stop`` statements, calls to ``mesa_error``, or
other error conditions.

A good test should run relatively quickly.  Costly parts can be
skipped over using a saved model when ``MESA_SKIP_OPTIONAL`` is set.

A good test should check its stopping condition before producing an
output model (i.e., set ``required_termination_code_string``).

A good test should use ``run_star_extras`` to perform more detailed
physical and/or numerical checks or report longitudinally interesting
values to the TestHub (see :ref:`developing/test_suite:MESA Test Hub`).

A good test should be numerically converged.  This is particularly
important in order to ensure that the test is robust.  Unconverged
tests can often fail in response to innocuous changes.  Models should
always be sufficiently converged that they run reliably and that any
quantities checked by the test are converged to better than the
tolerances of those checks.

.. note::

   In some instances, performing good science with a given test case
   may require even tighter convergence criteria that are practically
   excluded by test suite run time considerations.  In such a case,
   note this fact in the test case documentation.

A good test should have physically-motivated time step limits.  (This
often also an important part of ensuring convergence.)  The test
should trigger few retries and have few time steps limited by solver
convergence or by catch-all quantities like varcontrol.  The
``star_job`` options ``show_retry_counts_when_terminate = .true.`` and
``show_timestep_limit_counts_when_terminate = .true.`` will summarize
the retries and timestep limits encountered during the run.


Some tests should only be run in certain circumstances (e.g., if GYRE
is installed, if OP MONO opacities are present).  Such a test should
still compile and run even when these conditions are not met, but
should output the string "this test was intentionally skipped".  When
this string is present in the test output, the test scripts will
consider the test a success and no further checks will be performed.


README
~~~~~~

.. include:: ../../../star/test_suite/test_case_template/README.rst
    :start-line: 6



In the test suite directory
^^^^^^^^^^^^^^^^^^^^^^^^^^^

The existence of a sub-directory in ``test_suite`` is not sufficient
to tell MESA to perform a given test.  The list of tests associated with each module
lives in the file ``test_suite/do1_test_source``.  This has an a series of
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


In the docs directory
^^^^^^^^^^^^^^^^^^^^^

Once the ``README.rst`` file is created, a link in ``docs/source/test_suite`` should be created so that it will be rendered as part of the documentation.

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


Setting up new machine with MESA TestHub
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
1. Make an account on `the TestHub <https://testhub.mesastar.org/>`_ by messaging one of the admins (Bill Wolf). You will receive a token that you will need later in order to submit logs to the logs server. Remember also the email and password you use; you'll need these later in the terminal. 

2. Make sure ruby is installed, for example using the `Ruby version manager <https://rvm.io/>`_ (rvm).

  * If using the rvm, follow the instructions on that page to install the gpg keys. If this does not work, then execute the next line (``\curl -sSL https://get.rvm.io | bash -s stable``) instead, and follow the instructions printed in the terminal.
  
  * That line only installs the rvm; you also need Ruby itself. One can execute ``\curl -sSL https://get.rvm.io | bash -s stable --ruby`` as per the rvm installation page, but that requires sudo access. For a local installation, one can follow `this StackOverflow answer <https://stackoverflow.com/a/17219765>`_. 

3. Download and set up `mesa_test <https://github.com/MESAHub/mesa_test>`_, by doing the following:

  * ``gem install mesa_test``
  
  * ``mesa_test setup`` (here you will supply your email, password, and token from earlier; this will create a settings file in ``~/.mesa_test/config.yml``)
  
  * ``mesa_test install_and_test main`` will check out the main branch, test it, and submit the results to the testhub. 

4. If you want to set up ``mesa_test`` to run automatically on a cluster, Rob Farmer has created `a set of scripts <https://github.com/rjfarmer/mesa-helios-test>`_ that work with the Slurm workload manager. These scripts pull all commits and submit a job to the cluster queue for each new commit. You must edit the paths in all of the scripts to point to your own directories.

  * You can set up ``mesa_test`` or Rob's cluster script to run recurrently as a cronjob by doing ``crontab -e`` to edit the cronjob table. 
  
  * Add for example: ``10 * * * * ~/mesa/mesa-helios-test/runMesaTest.sh >/dev/null 2>&1`` to make it run every 10 minutes (or swap out ``runMesaTest`` with a ``mesa_test`` command). The parts at the end of that line prevent it from emailing you each time it runs.


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

Note that this string will also `prevent any GitHub actions from running <https://docs.github.com/en/actions/managing-workflow-runs/skipping-workflow-runs>`_.

[ci split]
^^^^^^^^^^

Splits the running of the test suite between machines. Current, if set, cannon will run the first half of the test cases while helios will run the second half.

[ci optional]
^^^^^^^^^^^^^

Runs ``MESA`` with the environment variable ``MESA_SKIP_OPTIONAL`` unset.  This requests that all parts of each test case be run (i.e., including optional parts).

[ci optional n]
^^^^^^^^^^^^^^^

Where ``n`` is an integer. Same as ``[ci optional]`` but only run the first ``n`` test cases.

[ci fpe]
^^^^^^^^

Compiles and runs ``MESA`` with the environment variable ``MESA_FPE_CHECKS_ON=1`` set. This requests that we turn on additional debugging checks.

[ci converge]
^^^^^^^^^^^^^

Runs the test suite with the environment variable ``MESA_TEST_SUITE_RESOLUTION_FACTOR`` set to a factor, giving a different temporal and spatial resolution (and max model number). 
