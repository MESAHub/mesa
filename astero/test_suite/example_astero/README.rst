.. _example_astero:

**************
example_astero
**************

An example optimisation run of the ``astero`` module, based on the
CoRoT target HD 49385.  This is the usual starting point if you want
to optimise model parameters using the ``astero`` module.  This test
case was introduced in Sec. 3.2 of |MESA II| and uses the
observational constraints presented by `Deheuvels et al. (2010)
<https://ui.adsabs.harvard.edu/abs/2010A%26A...515A..87D>`__ and
analysed by `Deheuvels & Michel (2011)
<https://ui.adsabs.harvard.edu/abs/2011A%26A...535A..91D>`__.  The
test uses the effective temperature, metallicity and mode frequencies
of angular degree 0, 1 and 2.

For speed, the test case only runs the initial parameters set by
``first_*``.  To optimse the model parameters, change ``search_type``
to one of the search options.

The test fails if the total |chi^2| is not between 1 and 15.  Mode
frequencies are quite sensitive to microphysics like opacity tables,
nuclear reaction rates and the equation of state, so changes to those
modules can eventually push |chi^2| above 15, triggering failure.
This is usually resolved by re-optimising the model.

If |chi^2| is smaller than zero, then |chi^2| was never successfully
evaluated.  The only possible negative value should be the placeholder
initial value of -1.

Finally, if |chi^2| is between zero and one, the test fails because
that sounds too good to be true!

Last-Update: 2020-11-13 (mesa r14889) by Warrick Ball
