.. _fast_simplex:

************
fast_simplex
************

Runs first few steps of a simplex optimisation for some synthetic data.
This test exists only to increase test coverage of the routines
involved in a simplex optimisation (some of which is used by other
optimisation routines).  It should not be the basis of any scientific
runs.  The tolerances are all very loose to make each iteration fast
enough that the test completes within a few minutes.  The data is taken from
the output of a 1 |Msun| model at an age around 1 Gyr.

The output in effect only checks that at least one set of trial
parameters returned a value of |chi^2|.  There are two main failure
modes.  First, failure might benignly indicate that the initial
parameters have strayed too far from the synthetic constraints to
produce any output, in which case the initial guesses or target data
should be adjusted.  Second, failure might indicate that the simplex
optimisation is genuinely broken and needs fixing.

Last-Update: 2022-02-06 (mesa 5805e64) by Warrick Ball
