.. _fast_from_file:

**************
fast_from_file
**************

Runs first few steps of the ``'from_file'`` search, which reads input
values from a data file and computes |chi^2| for each, based on some
synthetic data.  This test exists only to increase test coverage of
the routines.  It should not be the basis of any scientific runs.  The
tolerances are all very loose to make each iteration fast enough that
the test completes within a few minutes.  The data is taken from the
output of a 1 |Msun| model at an age around 1 Gyr.

To increase coverage, this test uses the custom parameters implemented
in ``my_param1``, ``my_param2`` and ``my_param3``.  Specifically, it
sets these to be the efficiency factors for semiconvection,
thermohaline mixing and gravitational settling.  Only ``my_param1``
and ``my_param3`` are varied.

The output in effect only checks that at least one set of trial
parameters returned a value of |chi^2|.  There are two main failure
modes.  First, failure might benignly indicate that the initial
parameters have strayed too far from the synthetic constraints to
produce any output, in which case the initial guesses or target data
should be adjusted.  Second, failure might indicate that the search
is genuinely broken and needs fixing.

Last-Update: 2022-03-29 (mesa 998d243) by Warrick Ball
