========================
Overview of rates module
========================

The ``rates`` module implements nuclear reaction rates.

Weak reactions
==============

The rates module provides MESA's weak reaction rates.

This includes a compilation of tables (in
``/data/rates_data/weakreactions.tables``) suitable for the high
densities and temperatures encountered in late stages of stellar
evolution.  These rates are based (in order of precedence) on the
tabulations of Langanke & Martı́nez-Pinedo (2000), Oda et al.  (1994),
and Fuller et al. (1985).  (These were referred to as a separate
``weaklib`` module in |MESA I|, but are now part of ``rates``.)

These reaction rates are tabulated as function of :math:`\rho Y_e` and
:math:`T`.  These tables cover 1 ≤ log(ρYe) ≤ 11 and 7 ≤ log(T) <
10.5, but are relatively coarse, with 11 points in the ρYe dimension
(∆ log ρYe = 1) and 12 points in the T dimension (∆ log T ≈ 0.25).
They are constructed assuming complete ionization (even if that is not
appropriate at the given conditions).

For each isotope pair, the tables have:

  * positron emission rate
  * electron capture rate
  * total neutrino energy loss rate
 
  * electron emission rate
  * positron capture rate
  * total anti-neutrino energy loss rate


These tables do not cover the lower density, lower temperature
conditions that are encountered in earlier stages of stellar evolution
(e.g., main sequence) or in stellar envelopes.  Their assumption of
complete ionization means that they are not physically appropriate for
conditions in which atoms are neutral or partially ionized (i.e., they
ignore electron captures from bound electrons).


We have a separate set of low temperature weak rates (internal
variable ``weak_lowT_rate``) that MESA will use when it is off the
weaklib tables.  The blend from these to the table occurs over a
specified temperature range (i.e., when ``T9_weaklib_full_on < T <
T9_weaklib_full_off``) and there is a separate blend for higher Z
elements (when ``Z >= weaklib_blend_hi_Z`` and
``T9_weaklib_full_on_hi_Z < T < T9_weaklib_full_off_hi_Z``).

.. warning::

  If you rely on any low temperature weak rates in your problem, you
  should carefully check them.  The current MESA approach often yields
  unphysical values, both in the low temperature/density limit and in
  the blending region.

Physically, the ``weak_lowT_rate`` should be the rate for a neutral
atom.  In practice (see ``rates/private/rates_initialize.f90``), these
rates are set by a call to ``reaclib`` at :math:`T = 10^7` K, or if
the rate is not in ``reaclib``, by consulting the file
``data/rates_data/weak_info.list``.  Depending on what information is
included in those files and on the physical processes that set the
decay, this approach may or may not yield the desired result.  (We
also must provide the average neutrino energies, which are not included in
``reaclib`` and the values in ``weak_info.list`` are mostly fantasies
or inappropriate.)




