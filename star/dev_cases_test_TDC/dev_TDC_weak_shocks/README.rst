Weak-shock pulse development case
=================================

This case updates the weak-shock pulse experiment distributed with the
MESA IV supporting inlists.  A localized energy source launches a pulse
through a 0.8 Msun white-dwarf envelope using ``u_flag`` hydrodynamics and
time-dependent convection.

The ``rn`` script runs the original five stages in order:

#. load and prepare the supplied starting model;
#. remove the central region;
#. settle the envelope with hydrodynamics enabled;
#. inject energy and follow the outgoing pulse; and
#. continue through the weak-shock breakouts.

The original numerical settings are retained where current MESA has a
direct equivalent.  Obsolete inlist controls and ``run_star_extras``
interfaces use their current forms.
