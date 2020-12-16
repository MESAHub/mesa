=====================
Environment Variables
=====================

The following lists environment variables that can effect the way MESA runs.

Essential
---------

MESA_DIR
~~~~~~~~

Location of the mesa source code.


GYRE_DIR
~~~~~~~~

Location of the gyre source code (only needed if running gyre).  If
you haven't moved it from the MESA_DIR then set as ``GYRE_DIR =
$MESA_DIR/gyre/gyre``.



Recommended
-----------

OMP_NUM_THREADS
~~~~~~~~~~~~~~~

Number of threads MESA will run with, should be between 1 and 2 * number of logical cores.


MESASDK_ROOT
~~~~~~~~~~~~

Location of the SDK.



Optional
--------

MESA_CACHES_DIR
~~~~~~~~~~~~~~~

Location where MESA will store and read cache files.
See star_job.defaults cache directories for more details.


MESA_TEMP_CACHES_DIR
~~~~~~~~~~~~~~~~~~~~

Location where mesa will write the cache file temporarily before moving to
MESA_CACHES_DIR. If set, this folder MUST be unique for each
MESA run. If not set defaults to ``./.mesa_temp_cache``.


MESA_OP_MONO_DATA_PATH
~~~~~~~~~~~~~~~~~~~~~~

Location of the OP_MONO data files
See $MESA_DIR/star/test_suite/radiative_levitation/inlist_radiative_levitation
for more details.


MESA_INLIST
~~~~~~~~~~~

By default MESA will look for a file called "inlist" in the local
working directory for its configuration. This overrides the filename
and can point to a file somewhere else.


Misc
----

Useful things that are not environment variables.

Command line arguments
~~~~~~~~~~~~~~~~~~~~~~

./star can accept one argument that overrides the location of the
inlist file (see mesa_inlist). This will also override the environment
variable mesa_inlist. ::

    ./star some_other_inlist_file

skip_build
~~~~~~~~~~

Empty file.  If present in $MESA_DIR then no compiling is done when ./install is invoked.
If present in a sub folder (e.g MESA_DIR/eos) then no compiling is done in the sub-folder

skip_test
~~~~~~~~~

Empty file.  Similar to skip_build, but this will skip the compile time tests.
If placed in a sub-folder then only the tests in that folder are skipped


Private
-------

These options are for developers, you do not need them for normal usage.
If you use them and things break then stop using them.

MESA_DIR_INTENTIONALLY_EMPTY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Bypass checks done at compile time for checking if MESA_DIR is set.


MESA_TEMP_CACHES_DISABLE
~~~~~~~~~~~~~~~~~~~~~~~~

If set, then we do not use the temp cache mechanism, writing cache
files directly to the $MESA_DIR/data/ folder.


MESA_ERROR_BACKTRACE_DISABLE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If set disables the generation of backtraces when we call mesa_error()
this is mostly helpful on macs as they don't generate useful backtraces.


MESA_FPE_CHECKS_ON
~~~~~~~~~~~~~~~~~~

When set to 1 this will turn on a series of compile time checks as well as
certain inlist options designed to catch floating point exceptions.
This should not be set during a normal run.
