Troubleshooting
===============

Getting help
------------


External packages (ADIPLS, GYRE, STELLA)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

MESA happily provides specific versions of external software
instruments such as GYRE, ADIPLS, and STELLA.  While the MESA
developers do not actively support these external packages, we do
support interfaces to them.  See, for example, 
``astero_adipls`` and ``astero_gyre`` within ``$MESA_DIR/astero/test_suite`` 
for how to make calls to those packages
during a stellar evolution run. Also see
``$MESA_DIR/star/test_suite/ccsn_IIp`` for how to prepare input for STELLA
from the output of mesa/star run. Finally,
``$MESA_DIR/stella/res/stella_extras.f90`` has routines for converting some of
the output of STELLA into a more "MESA friendly" form. Again, we
support all of these interfaces and continue to welcome questions
about using them and suggestions for making them better.

If the above instruments are useful for one's experiments, then
fantastic. If one encounters difficulties or bugs when using one of
these packages, then (a) don't use the package in the troublesome
regime and (b) direct reports to the developers of those
packages. GYRE is an excellent example. One of reasons why one doesn't
see many GYRE related questions on mesa-users is because GYRE has an
active and robust support system (see
http://www.astro.wisc.edu/~townsend/gyre-forums/). This said, it is
possible that another mesa-user has encountered the same challenge
when using one of these external instruments, been able to resolve the
issue, and is willing to share their solution with the community. We
continue to welcome mesa-users posts about these external packages.
