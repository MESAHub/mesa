=====================
Environment Variables
=====================

The following lists environment variables that can effect the way MESA runs.

Essential
---------

MESA_DIR
~~~~~~~~

Location of the MESA source code.


GYRE_DIR
~~~~~~~~

Location of the GYRE source code (only needed if running GYRE).  If
you haven't moved it from the ``$MESA_DIR`` then set as 
``GYRE_DIR=$MESA_DIR/gyre/gyre``.



Recommended
-----------

OMP_NUM_THREADS
~~~~~~~~~~~~~~~

Number of threads MESA will run with. Should be between 1 and 2 * number of logical cores.


MESASDK_ROOT
~~~~~~~~~~~~

Location of the SDK.



Optional
--------

MESA_CACHES_DIR
~~~~~~~~~~~~~~~

Location where MESA will store and read cache files.
See :ref:`star_job.defaults cache directories 
<reference/star_job:cache directories>` for more details.
Beware when using this and changing MESA versions; 
best practice would be to remove its contents upon changing versions. 


MESA_TEMP_CACHES_DIR
~~~~~~~~~~~~~~~~~~~~

Location where MESA will write the cache file temporarily before moving to
``$MESA_CACHES_DIR``. If set, this folder MUST be unique for each
MESA run. If not set defaults to ``./.mesa_temp_cache``.


MESA_OP_MONO_DATA_PATH
~~~~~~~~~~~~~~~~~~~~~~

Location of the OP_MONO data files.
See the :ref:`radiative levitation <radiative_levitation>` test suite 
(located in ``$MESA_DIR/star/test_suite/radiative_levitation/inlist_radiative_levitation``)
for more details.


MESA_INLIST
~~~~~~~~~~~

By default MESA will look for a file called ``inlist`` in the local
working directory for its configuration. This overrides the filename
and can point to a file somewhere else.


MESA_FORCE_PGSTAR_FLAG
~~~~~~~~~~~~~~~~~~~~~~

If set, this shell variable will override the inlist variable ``pgstar_flag``.

If set to either ``true`` or ``TRUE`` MESA will act as if the user had set ``pgstar_flag=.true.`` in their inlist.
If set to either ``false`` or ``FALSE`` MESA will act as if the user had set ``pgstar_flag=.false.`` in their inlist.
If set to anything else (or not set at all) MESA will use the value as set in the inlists.


MESA_SKIP_OPTIONAL
~~~~~~~~~~~~~~~~~~

If set then when running a test_suite case, then certain optional inlists will be skipped. 
If not set, then all inlists will be run.



Misc
----

Useful things that are not environment variables.

Command line arguments
~~~~~~~~~~~~~~~~~~~~~~

``./star`` can accept one argument that overrides the location of the
inlist file. This will also override the environment variable ``$MESA_INLIST``. ::

    ./star some_other_inlist_file

skip_build
~~~~~~~~~~

Empty file.  If present in ``$MESA_DIR`` then no compiling is done when ``./install`` is invoked.
If present in a sub folder (e.g ``$MESA_DIR/eos``) then no compiling is done in the sub-folder.

skip_test
~~~~~~~~~

Empty file.  Similar to ``skip_build``, but this will skip the compile time tests.
If placed in a sub-folder, then only the tests in that folder are skipped. 


Private
-------

These options are for developers; you do not need them for normal usage.
If you use them and things break then stop using them.

MESA_DIR_INTENTIONALLY_EMPTY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Bypass checks done at compile time for checking if ``$MESA_DIR`` is set.


MESA_TEMP_CACHES_DISABLE
~~~~~~~~~~~~~~~~~~~~~~~~

If set, then we do not use the temp cache mechanism, writing cache
files directly to the ``$MESA_DIR/data/`` folder.


MESA_ERROR_BACKTRACE_DISABLE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If set, disables the generation of backtraces when we call ``mesa_error()``.
This is mostly helpful on macs, as they don't generate useful backtraces.


MESA_FPE_CHECKS_ON
~~~~~~~~~~~~~~~~~~

When set to 1 this will turn on a series of compile time checks as well as
certain inlist options designed to catch floating point exceptions.
This should not be set during a normal run.


MESA_TEST_SUITE_RESOLUTION_FACTOR
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If set to a value other than 1, then ``mesh_delta_coeff`` and 
``time_delta_coeff`` will be multiplied by its value, and 
``max_model_number`` will be multiplied by its inverse. 
For example, this can be set to 0.5 to double the space and time resolution 
as well as the maximum model number. 
During a normal run, this should be either set to 1, or not set. 
