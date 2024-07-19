*********
Changelog
*********

.. warning:: As of r24.08.1, building MESA now requires Python (3.5 or newer) be installed.

Changes in main
===============

.. note:: This describes changes present in the development version of MESA (``main`` branch) relative to the most recent release.

.. _New Features main:

New Features
------------

**Replacement of HDF5io by ForUM**

The HDF5io module (in the :file:`hdf5io` subdirectory) has been
replaced by the ForUM module (in the :file:`forum`
subdirectory). HDF5io provided a high-level wrapper around the `HDF5
input/output library <https://www.hdfgroup.org/solutions/hdf5/>`__; it
was based on source files copied from an alpha version of
ForUM. :git:`ForUM <rhdtownsend/forum>` (short for Fortran Utility
Modules) is a small library providing input/output, operating system,
memory management and utility routines for Fortran.

Now that ForUM has officially been released, and also that recent GYRE
releases also require ForUM, it made sense to include the official
distribution of ForUM within MESA. One important consequence of this
decision is that ForUM uses :git:`fypp <aradi/fypp>` for
pre-processing and templating, and fypp in turn requires a Python
interpreter. As a consequence, building MESA now requires Python (3.5
or newer) be installed.

**Update to GYRE 7.2.1**

The GYRE distribution bundled with MESA has been updated to release
7.2.1. Full details about this release can be found on the `GYRE
documentation page <https://gyre.readthedocs.io/en/v7.2.1/>`__. From
the perspective of MESA users, the significant change arising from
this update is that the MESA-to-GYRE interface is now provided via the
library file :file:`libgyre_mesa.a` and the module file
:file:`gyre_mesa_m.mod` (previously, these were :file:`libgyre.a` and
:file:`gyre_lib.mod`, respectively). A new variable, `LOAD_GYRE`, has
been added to :file:`utils/makefile_header` to simplify linking
against :file:`libgyre_mea.a`. These changes will likely only affect
those users that make calls to GYRE from inside
:file:`run_star_extras.f90`.

.. _Bug Fixes main:

Bug Fixes
---------


Changes in r24.06.1-rc1
=======================

.. _New Features r24.06.1-rc1:

New Features
------------

``max_allowed_nz`` is now ignored if the value is less than or equal to zero.

Kap
~~~

**High Temperature Opacity Tables**

Type 1 Rosseland-mean opacity tables from The Los Alamos
OPLIB database (`Colgan et al. 2016 <https://ui.adsabs.harvard.edu/abs/2016ApJ...817..116C/abstract>`_) are now available (`Farag et al. 2024 <https://doi.org/10.3847/1538-4357/ad4355>`_).
These tables cover the region :math:`0.0 \leq X \leq 1-Z` and
:math:`0.0\leq Z \leq 0.2`. Each set of OPLIB
opacity tables contains 1194 individual tables, a
dramatic increase in table density (in the X--Z plane) over the standard
126 individual tables provided in previous opacity releases. These tables
are available for four solar-scaled abundance mixtures constructed from photospheric estimates of the solar heavy element abundance:
(`GS98, Grevesse & Sauval 1998 <https://ui.adsabs.harvard.edu/abs/1998SSRv...85..161G/abstract>`_),
(`AGSS09(a09p), Asplund et al. 2009 <https://ui.adsabs.harvard.edu/abs/2009ARA%26A..47..481A/abstract>`_),  
(`AAG21, Asplund et al. 2021 <https://ui.adsabs.harvard.edu/abs/2021A%26A...653A.141A/abstract>`_),
and (`MB22, Magg et al. 2022 <https://doi.org/10.1051/0004-6361/202142971>`_). Users can
adopt this new set of tables by selecting one of the following
options for ``kap_file_prefix``:

+ ``'oplib_gs98'``
+ ``'oplib_agss09'``
+ ``'oplib_aag21'``
+ ``'oplib_mb22'``
 
See :ref:`kap/overview:Overview of kap module` and
:ref:`kap/defaults:kap_file_prefix` for more details on the
implementation of these tables. For further details on these new OPLIB opacity tables, a direct comparison with 
the Type 1 OPAL/OP tables as well as their effect on solar models can be found in
in `Farag et al. 2024 <https://doi.org/10.3847/1538-4357/ad4355>`_.


**Low Temperature Opacity Tables**

Low temperature Rosseland-mean opacity tables for both (`AAG21, Asplund et al. 2021 <https://ui.adsabs.harvard.edu/abs/2021A%26A...653A.141A/abstract>`_),
and (`MB22, Magg et al. 2022 <https://doi.org/10.1051/0004-6361/202142971>`_) 
solar-scaled abundance mixtures have been privately communicated by Jason Ferguson. These opacity tables were
computed following the approach of `Ferguson et al. (2005) <https://ui.adsabs.harvard.edu/abs/2005ApJ...623..585F/abstract>`_. Users can
adopt this new set of tables by selecting one of the following
options for :ref:`kap/defaults:kap_lowT_prefix`:

+ ``'lowT_fa05_mb22'``
+ ``'lowT_fa05_aag21'``

**Opacity interpolation**

We have updated the opacity interpolation scheme to provide much higher quality derivatives when doing cubic interpolation in composition.

MESA interpolates across opacity tables in the :math:`X–Z` plane through the use of two consecutive 1D splines.
MESA offers users the ability to choose linear or cubic interpolation for these splines, 
while leaving the default as linear interpolation::

  cubic_interpolation_in_X = .false.
  cubic_interpolation_in_Z = .false.

This choice of default was primarily due to the fact that
the previous cubic composition interpolation scheme in MESA suffered from poor quality interpolated opacity derivatives with respect to
density and temperature, which often disagreed with the numerical derivatives produced via nearest neighbor
Richardson extrapolation. The figure below shows this comparison on a logarithmic scale, where in general red indicates poor quality
derivatives and blue indicates high quality derivatives.

.. figure:: changelog_plots/cubic_dfridr_dkapdT.png
   :alt: old cubic relative kap derivative error

   This figure shows the logarithmic relative error in the derivative :math:`\partial \kappa / \partial T` (:math:`X` = 0.625, :math:`Z` = 0.015),
   for an OPAL opacity table grid using Grevesse & Sauval (1998) abundances, generated from MESA’s kap module, using the previous cubic interpolation scheme.
   The OPLIB log(:math:`R`) = −8, 1.5 table boundaries are marked with a solid black line and the OPAL/OP log(:math:`R`) = 1.0 boundary is shown with a dashed line.
   The approximate location of the Z-dependent transition to an electron conduction dominated opacity is marked with dot-dash blue curve.
   Regions for Atomic, molecular, and compton scattering opacity are labeled and presented with their associated blending regions.


While the opacity derivatives do not directly appear in the canonical equations of stellar structure, they do appear in the Jacobian matrix for MESA's implicit solver.
Numerically unstable opacity derivatives can halt the progress of the solver and ultimately crash a calculation. 

To improve the numerical stability of MESA's cubic opacity interpolation routines, we have implemented
automatic differentiation into the opacity interpolating functions. Now, when using cubic interpolation, the opacity derivatives for an arbitrary mixture
in the :math:`X–Z` plane are computed by taking the derivative of the interpolating function as opposed to the interpolant of the derivatives. This improvement
has led to a significant reduction in the relative derivative error and an increase in the numerical accuracy of opacity derivatives computed with cubic interpolation. 

.. figure:: changelog_plots/cubic_dfridr_dkapdT_ad.png
   :alt: new cubic relative kap derivative error

   Same as previous figure, but for new cubic interpolation scheme taking advantage of automatic differentiation.


This new implementation of cubic interpolation in composition for opacity tables comes close to achieving the derivative quality of the linear interpolation
option (shown below), while also providing more accurate opacity physics between opacity table grid points.

.. figure:: changelog_plots/linear_dfridr_dkapdT.png
   :alt: linear relative kap derivative error

   Same as previous figure, but for linear interpolation instead of cubic.


For this MESA release, linear interpolation remains the default method for interpolating in composition between opacity tables
while we continue to investigate the residual areas where cubic interpolation appears to occasionally produce lower quality derivatives.
However, adopting cubic interpolation has been shown to consistently increase the overall 
opacity of a model, and can directly effect the structure of solar models, see Appendix B & C in `Farag et al. 2024 <https://doi.org/10.3847/1538-4357/ad4355>`_.
We anticipate making cubic interpolation the default in a future MESA release version. 
We encourage users to experiment with these different opacity interpolation routines and be mindful of the effect they can have on their stellar models.


Chem
~~~~
New initial metal mass fractions ``initial_zfracs`` taken from photospheric estimates of the solar heavy element abundances in (AAG21, Asplund et al. 2021) and (MB22, Magg et al. 2022)
are now available. See :ref:`reference/star_job:initial_zfracs` for more details.

.. _Bug Fixes r24.06.1-rc1:

Bug Fixes
---------

Rates
~~~~~

There has been a bug present in the rates module due to the incorrect
phase space factors for reverse reaction rates involving greater than 2 reactants or 
products. This bug resulted in inconsistent equilibrium compositions when the network
evolves into nuclear statistical equilibrium (NSE), at temperatures exceeding 4 GK. 
This bug effects users who evolve models into NSE using large reaction networks. This
includes evolving massive stars to core-collapse. Smaller networks such as the ``approx21``
networks are less affected. We strongly recommend that users update to the latest MESA release.

See `gh-575 <https://github.com/MESAHub/mesa/issues/575>`_

``max_allowed_nz``
~~~~~~~~~~~~~~~~~~

MESA no longer produces a segmentation fault if it tries to increase
the number of cells beyond ``max_allowed_nz``.

``pgbinary``
~~~~~~~~~~~~
A bug in the `pgbinary` axes definitions resulted in models crashing when running in single star mode and has been corrected.
Another bug inhibited the mass of the not modeled star from being displayed in the `pgbinary` panel and has also been corrected.

See `gh-634 <https://github.com/MESAHub/mesa/issues/634>`_ 


Changes in r24.03.1
===================

.. _Backwards-incompatible changes r24.03.1:

Backwards-incompatible changes
------------------------------

Hooks
~~~~~

The ``other_pressure`` hook has been converted to use ``auto_diff``
thus the variable ``s% extra_pressure`` is now an ``auto_diff``
and allows for the setting of the partial derivatives the
pressure with respect to other variables.

run_star_extras
~~~~~~~~~~~~~~~

Previously we had logic to determine if an extra history value should be saved
as an int or a float (users can only provide data as a float). This was error
prone, so now we save extra history values as floats.

Asteroseismic scaling relations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The default values of the parameters ``delta_nu_sun`` and ``nu_max_sun`` used in the asteroseismic scaling relations
have been updated to the values reported by `Lund et al. (2017) <https://ui.adsabs.harvard.edu/abs/2017ApJ...835..172L/abstract>`_
(:math:`\Delta\nu_\odot = \mathrm{134.91}\,\mu\mathrm{Hz}` and :math:`\nu_{\mathrm{max},\odot} = \mathrm{3078}\,\mu\mathrm{Hz}`).
The parameter ``Teff_sun`` has been renamed ``astero_Teff_sun`` to avoid confusion
with the fixed constant ``Teffsun`` in the ``const`` module.  ``astero_Teff_sun``'s
default value has been set to 5772 K, as in `IAU 2015 Resolution B3 <https://ui.adsabs.harvard.edu/abs/2015arXiv151007674M>`_.

stash.py
~~~~~~~~
The script ``stash.py`` for storing MESA runs has been moved out of ``$MESA_DIR/star/work``
and into the ``$MESA_DIR/scripts``. Consider adding this directory to your path. 


.. _New Features r24.03.1:

New Features
------------

shmesa
~~~~~~

We have introduced a new set of command line utilities for interacting with MESA. 
See the README in ``$MESA_DIR/scripts/shmesa``, or online `here <https://github.com/MESAHub/mesa/tree/main/scripts/shmesa>`_. 

These utilities provide functionality such as changing inlist parameters (``shmesa change``) or filling in the full 
``run_star_extras.f90`` template (``shmesa extras``). 

We recommend adding ``shmesa`` to your ``PATH`` (via, e.g., placing ``PATH=$PATH:$MESA_DIR/scripts/shmesa`` in your ``~/.bash_profile``). 


Rates
~~~~~

In ``$MESA_DIR/data/rates_data/rate_tables`` we now ship thirteen C12(a,g)O16
rates spanning the uncertainty range of -3 to +3 sigma from `Deboer et al. 2017 <https://ui.adsabs.harvard.edu/abs/2017RvMP...89c5007D/abstract>`_
with updated numerical resolution from `Mehta et al. 2022 <https://ui.adsabs.harvard.edu/abs/2022ApJ...924...39M/abstract>`_. 

These rates can be accessed through the rate selection mechanism. 
These can be loaded, either by the normal mechanism of adding the filename
to a ``rates_list`` file, or by using the option ``filename_of_special_rate``.
Several examples in the test suite now make use of these rates, such as
massive stars and models for building white dwarfs.

Maximum net size
~~~~~~~~~~~~~~~~

Previous MESA versions where limited to ~300 isotopes in the nuclear network,
this has now been (slightly) relaxed.
In principle the maximum network size is sqrt(2^31/nz). Thus for nz=1000 zones
this is about ~1400 isotopes, at 5000
zones this is ~650 isotopes. Note this will require a huge amount of RAM,
``mesa_495.net`` requires at least ~80GB of RAM, and it is estimated
that for 1400 isotopes you will need ~650GB of RAM.

Björklund Wind
~~~~~~~~~~~~~~

The Björklund et al. (2021) wind scheme has been implemented for use in all wind
schemes. E.g. ``hot_wind_scheme = 'Bjorklund'`` in ``&controls``.

GYRE Update
~~~~~~~~~~~

The bundled GYRE oscillation code has been updated to `release 7.1
<https://gyre.readthedocs.io/en/v7.1/>`_. The GYRE and GSM pulsation
output formats have been also been updated to take advantage of new
features in GYRE. The version of an output file can be set by the new
``gyre_data_schema`` parameter in ``&controls``; the default value is 101 (version 1.01)
produced GYRE-format files that are backward-compatible with older
GYRE releases. (If using GSM-format files, set to 110 instead for
backward compatibility).

White Dwarf O/Ne Phase Separation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Phase separation for crystallizing white dwarfs now supports options
for either C/O or O/Ne phase diagrams. You can select the appropriate
option with the ``phase_separation_option`` in the ``&controls`` inlist
section (see documentation at :ref:`reference/controls:phase_separation_option`).
The phase diagram for O/Ne separation comes from
`Blouin & Daligault (2021b) <https://ui.adsabs.harvard.edu/abs/2021ApJ...919...87B/abstract>`_.

Massive Star test_suite Updates
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ``12Msun_pre_ms_to_core_collapse``, ``20Msun_pre_ms_to_core_collapse``, and ``zams_to_cc_80``
test_suites have each been updated and now fully evolve models at solar metallicities to core
collapse while keeping their surface boundaries at tau = 1. The models offer a framework which
has been tested to function with large reaction networks containing > 200 isotopes thanks to
'op_split_burn'. This updates includes the revival of the ``make_pre_ccsn_13bvn`` test_suite
as seen in section 6.9 of MESA IV.

Reintroduction of Velocity Drag for v_flag
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The option to apply a drag term to the velocities in the outer envelope
with the control ``drag_coefficient`` and ``min_q_for_drag`` was removed after
version 15140, but has now been reintroduced. The logical control ``use_drag_energy``
has been introduced as well, so users can decide whether the drag energy is included
in the energy equation. We also provide an additional option to turn off the 
velocities above a desired optical depth with the control ``velocity_tau_lower_bound``.

.. _Bug Fixes r24.03.1:

Bug Fixes
---------

The ZAMS model data used to generate the initial model for simulations with
``create_pre_main_sequence_model = .false.`` and ``load_saved_model = .false.``
had an issue where stars between 1.0 and 1.58 Msun would have a starting 
central hydrogen mass fraction markedly lower than other masses. The grid of 
starting models has been recalculated with a more stringent stopping condition, 
and now all pre-computed ZAMS models have a central hydrogen mass fraction very 
near 0.697.

The ``fixed_Teff``, ``fixed_Tsurf``, ``fixed_Psurf``,  and ``fixed_Psurf_and_Tsurf``
atmosphere options were removed in r15140. We have reimplemented them although we
caution users that their implementation could conflict with ``mlt_option = 'TDC'``.

The EOS coverage regions have been updated to fall back to ideal gas in a region
previously covered by HELM where it returned unphysical floor values of ``1e-20``
for pressure, energy, and entropy. The most up-to-date EOS coverage plots can
be found in the EOS documentation: :ref:`eos/overview:Overview of eos module`.

A bug in the implementation of the ``fe_core_infall_limit`` sometimes resulted 
in the premature termination of a model due to large negative velocities outside
the fe core.

 See `gh-626 <https://github.com/MESAHub/mesa/issues/626#issue-2157840823>`_

Changes in r23.05.1
===================

.. _Backwards-incompatible changes r23.05.1:

Backwards-incompatible changes
------------------------------

Extra inlist controls are now arrays
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Almost all MESA inlists have the option of reading other inlists,
which is a feature canonically used in the main ``inlist`` file.
e.g. the ``inlist`` in the standard ``star/work`` directory has ::

    read_extra_controls_inlist1 = .true.
    extra_controls_inlist1_name = 'inlist_project'

where the inlist number could range from 1 to 5.

These and all similar controls have been replaced with arrays like ::

    read_extra_controls_inlist(1) = .true.
    extra_controls_inlist_name(1) = 'inlist_project'

That is, the number should be moved to the end of the control name
and placed in round brackets.

This allows a lot of duplicate code to be refactored but will break
almost all existing MESA inlists.  To update an old inlist to this new
style, you can use the following ``sed`` terminal command: ::

    sed -E '/inlist[1-5]/s/([1-5])([_a-z]*) =/\2\(\1\) =/' -i inlist_name

where ``inlist_name`` is the inlist (or inlists) that you'd like to update.
This will *replace* the file ``inlist_name``.  Omit the ``-i`` flag if you'd
like to see the changes without modifying the file.

``sed`` is a standard tool that is included with macOS and most Linux distributions.
For convenience, we have also included a bash script that will call a version of
this ``sed`` command (along with ``sed`` commands for the next changlog entry as well)
to update all inlist files (``inlist*``), which you can run in any work directory
where you want to update every inlist by invoking ::

  $MESA_DIR/scripts/update_inlists

This script will save the previous versions of your inlists to a directory named
``backup_inlists``.
  
Renamed controls for upper limits
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The following controls in ``&controls`` for upper limits on
`when to stop <reference/controls.html#when-to-stop>`_ have
been renamed:

+------------------------------+------------------------------------+
+ Old                          + New                                +
+==============================+====================================+
+                              +                                    +
+ ``log_center_density_limit`` + ``log_center_density_upper_limit`` +
+                              +                                    +
+ ``log_center_temp_limit``    + ``log_center_temp_upper_limit``    +
+                              +                                    +
+ ``center_entropy_limit``     + ``center_entropy_upper_limit``     +
+                              +                                    +
+ ``max_entropy_limit``        + ``max_entropy_upper_limit``        +
+                              +                                    +
+------------------------------+------------------------------------+

You can substitute the new names for the old ones using the command
line tool ``sed`` with, e.g. ::

    $ sed 's/log_center_density_limit/log_center_density_upper_limit/' -i <inlist_filename>

Abundance-based timestep controls are now arrays
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The previous controls ::

    dH_limit_min_H = 1d99
    dH_limit = 1d99
    dH_hard_limit = 1d99
    dH_decreases_only = .true.
    dH_div_H_limit_min_H = 1d-3
    dH_div_H_limit = 0.9d0
    dH_div_H_hard_limit = 1d99

and similar controls for ``He`` and ``He3`` have been replaced with
arrays.  This simplifies the code and allows the controls to be
applied to any species in the net.  A new control
``dX_limit_species(...)`` specifies which elements will be checked.
The previous behaviour for ``H``, for example, has been replaced with
::

    dX_limit_species(1) = 'h1'
    dX_limit_min_X(1) = 1d99
    dX_limit(1) = 1d99
    dX_hard_limit(1) = 1d99
    dX_decreases_only(1) = .true.
    dX_div_X_limit_min_X(1) = 1d-3
    dX_div_X_limit(1) = 0.9d0
    dX_div_X_hard_limit(1) = 1d99

The old ``H``, ``He`` and ``He3`` controls correspond to species
``h1``, ``he4`` and ``he3``, respectively.  You can also set the
species to ``X``, ``Y`` or ``Z``, in which case the checks are applied
*individually* to all isotopes of hydrogen, helium or metals,
respectively.

Colors
~~~~~~

The colors module now returns ``-1d99`` when asking for a value that is
off table.


.. _New Features r23.05.1:

New Features
------------

White Dwarf C/O Phase Separation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

An option to include carbon-oxygen phase separation for crystallizing C/O white dwarfs is now available,
using the phase diagram of `Blouin et al. (2021) <https://ui.adsabs.harvard.edu/abs/2021PhRvE.103d3204B/abstract>`_.
The MESA implementation is described in `Bauer (2023) <https://ui.adsabs.harvard.edu/abs/2023arXiv230310110B/abstract>`_.
More documentation and associated controls can be found at :ref:`reference/controls:do_phase_separation`.
This option is off by default, but it is on in the ``wd_cool_0.6M`` test case.

Module enhancement: ``pgbinary``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

When running ``./binary`` models it is useful to have graphical
output to 'see' what's going on.
Previously, this was only possible on the ``pgstar`` level, meaning you would
need to setup two ``pgstar`` windows if you are evolving two stars in the
binary.

Here we introduce ``pgbinary``, which acts much like ``pgstar``. You enable it
with the ``&binary_job`` inlist with ``pgbinary_flag = .true.``. Then you
select windows and/or files to be plotted in the ``&pgbinary`` inlist.
Currently the following plot types can be created:

* History_Track[1-9],
* Summary_History,
* History_Panels[1-9],
* Text_Summary[1-9],
* Grid[1-9],

analogous to their ``pgstar`` equivalents, and two ``pgbinary``-only plots:

* Star[1-2], to plot a star window through ``&pgstar`` controls, within ``pgbinary``.
* Orbit, a visual representation of the stars' sizes to their separation

Main use case is to have a single window containing both stars' ``pgstar`` info,
through using ``Grid`` at the ``pgbinary`` level, populating it with ``Star1``
and ``Star2``, and have each plot profile info, Kipp diagrams etc...

Resolution control convective_bdy_weight has been reintroduced
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The option to add extra resolution at convective boundaries
with the control ``convective_bdy_weight`` was removed after
version 12115, but has now been reintroduced in a simplified form.
This control no longer applies to newly nonconvective zones, but
does add resolution at the location of convective boundaries.
This was found to be important for smooth convective boundary
evolution with convective premixing.

Hooks
~~~~~

A new other_close_gaps hook has been added. Provided by Simon Guichandut

.. _Bug Fixes r23.05.1:

Bug Fixes
---------

Rates
~~~~~

There has been a bug present in the rate ``r_c12_to_he4_he4_he4`` in r22.05.1 and r22.11.1. 
This causes an excessive amount of C12 to be burnt during core helium burning. 
We strongly recommend that users update to the latest MESA.

See `gh-526 <https://github.com/MESAHub/mesa/issues/526>`_

RTI
~~~

A bug has existed since shortly after r15140 where RTI
mixing will be effectively zero in a model even with the
``RTI_flag=.true.`` set.

This has now been fixed. Users of RTI mixing are recommended to upgrade to the
newest MESA version.

See `gh-503 <https://github.com/MESAHub/mesa/issues/503>`_


Changes in r22.11.1
===================

.. _Backwards-incompatible changes r22.11.1:

Backwards-incompatible changes
------------------------------

.. note::

   A large amount of internal clean up has occurred since the last release.  This lists some of the most important changes, but the list is not exhaustive.


Module-level changes
--------------------

astero
------

The main controls for the selection of parameters and non-seismic
constraints (which were dubbed "variables") has changed.  The defaults
files document the new interface but the most important changes are
repeated here.

Each non-seismic constraint is now given a name, target value,
uncertainty and flag for whether to include it in the |chi^2|
calculation.  The default ``work`` folder will either try for one of
the custom options included for backward compatibility (e.g. ``Rcz``)
or fall back to computing the matching history column (e.g. for
``log_g``).  So whereas an effective temperature constraint would
previously be included using, say ::

    include_Teff_in_chi2_spectro = .true.
    Teff_target = 6000
    Teff_sigma = 100

you would now use ::

    constraint_name(1) = 'Teff'
    include_constraint_in_chi2_spectro(1) = .true.
    constraint_target(1) = 6000
    constraint_sigma(1) = 100

The maximum number of such constraints is currently 100 but can
trivially be increased at compile time by modifying
``max_constraints`` in ``astero/public/astero_def.f90``.

Similarly, each parameter now has a name, initial value, minimum,
maximum and grid-spacing.  So whereas the mixing-length parameter
was previously controlled with something like ::

    vary_alpha = .true.
    first_alpha = 1.7
    min_alpha = 1.5
    max_alpha = 1.9
    delta_alpha = 0.1

you would now use ::

    param_name(1) = 'alpha'
    vary_param(1) = .true.
    first_param(1) = 1.7
    min_param(1) = 1.5
    max_param(1) = 1.9
    delta_param(1) = 0.1

Again, the maximum number of parameters is 100 and can be increased at
compile time by modifying ``max_parameters`` in
``astero/public/astero_def.f90``.

The default ``run_star_extras.f90`` defines the hooks
``set_constraint_value`` and ``set_param`` so that the old options
remain available, though with a new syntax.  Users can also use those
routines to define their own parameters and constraints.

The output files contain information for constraints or parameters
with names that are not ``''``.  Thus, the column order now varies but
the same information is present and now follows the same structure as
histories and profiles.

star
----

For wind mass loss, schemes that scale with metallicity now employ
``Zbase`` rather than ``Z`` (as long as ``Zbase`` is set to a non-negative value,
otherwise we fall back to ``Z``). This reflects the fact that wind recipes primarily
account for the opacity of iron-group elements, which have surface abundances that
are unlikely to change during evolution. This change therefore avoids
unphysical influence on winds by, e.g., evolution of surface CNO abundances.

test_suite
----------

All test cases have now had the inlist option, makefile variable, and shell script variable, ``MESA_DIR`` removed.
This means that you no longer need to do anything to make a ``MESA`` test case work outside of the test suite.
Test cases now infer their ``MESA_DIR`` variable entirely by the environment variable ``$MESA_DIR``. 

The history output option ``tri_alfa`` (and other quantities that relate to the triple-alpha nuclear reaction) have been renamed to ``tri_alpha`` for better consistency with other ``_alpha`` reactions.

net
---

The derived type net_info (conventional given the symbol ``n``) is no longer a pointer. If you declare a local copy of the variable, you should also ensure to do ``n% g => g`` to make sure that net_info knows
about the ``net_general_info`` derived type. ``g`` can be had from a call to ``get_net_ptr(handle, g, ierr)``.

The pointer array ``net_work`` and its size ``net_lwork`` have been removed from the net interface, thus these variables should be removed form any ``other_net_get`` and ``other_split_burn`` hooks.
The following routines have also been removed as they are no longer needed ``net_work_size``, ``get_net_rate_ptrs``, ``net_1_zone_burn_work_size``, ``get_burn_work_array_pointers``, ``net_1_zone_burn_const_density_work_size``, and ``get_burn_const_density_work_array_pointers``

Previously you could pass ``arg_not_provided`` for either the temperature (density) or log(temperature) (log(density)). Now you must pass both explicity.

ADIPLS
------

ADIPLS now has a ``USE_ADIPLS`` flag in ``utils/makefile_header`` to enable is build to be disabled.


Changes in r22.05.1
===================

.. _Backwards-incompatible changes r22.05.1:

Backwards-incompatible changes
------------------------------

.. note::

   A large amount of internal clean up has occurred since the last release.  This lists some of the most important changes, but the list is not exhaustive.


Module-level changes
--------------------

astero
------

``&astero_search_controls`` now has an option ``astero_results_directory`` to specify a folder into which all of
``astero``'s results are saved (like ``log_directory`` in ``star``).  The default is ``outputs``, so if you
can't seem to find your ``astero`` output, have a look there.

``&astero_search_controls`` also now has options ::

    astero_results_dbl_format = '(1pes26.16)'
    astero_results_int_format = '(i26)'
    astero_results_txt_format = '(a26)'

by which the user can set the formats of floats, integers and strings in the ``astero`` results file,
much like ``star_history_*_format`` does for history files.

The format of the ``astero`` results file has changed to match histories and profiles.
The contents of the file are unchanged.

rates
-----

The 7Be(e-,nu)7Li has been switched from REACLIB rate to that of `Simonucci et al 2013 <https://ui.adsabs.harvard.edu/abs/2013ApJ...764..118S/abstract>`_. This is
due to the fact that the REACLIB rate does not take into account the neutral ion rate below 10^7 K.

The ability to set the rates preferences has been removed. This added alot of complexity to the rates code handling NACRE and REACLIB and made it difficult to reason about where a rate actually came from.
From now on we excusivily use NACRE for any rate that cares about temperatures below 10^7K (for all temperatures), REACLIB for almost all other rates, and a small number of rates
from CF88 (if they aren't in REACLIB or NACRE). 

Of note is that the default C12(a,g)O16 rate has thus changed from NACRE to that of REACLIB.

The options ``set_rates_preferences``, ``new_rates_preference``, and ``set_rate_c1212`` have been removed without replacements.

The options ``set_rate_c12ag``, ``set_rate_n14pg``, and ``set_rate_3a`` have been removed. However, those rates can now be access thorugh a new 
rate selection mechanism. In ``$MESA_DIR/data/rates_data/rate_tables`` we now ship a number of MESA provided rate tables. These can be loaded,
either by the normal mechanism of adding the filename to a ``rates_list`` file, or by using the new option ``filename_of_special_rate``.
This option sets the filename to load the rate from for the rate specified by ``reaction_for_special_factor``.

Thus the options:

::

      num_special_rate_factors = 1
      reaction_for_special_factor(1) = 'r_c12_ag_o16'
      special_rate_factor(1) = 1
      filename_of_special_rate(1) = 'r_c12_ag_o16_kunz.txt'

replaces the previous:

::

  set_rate_c12ag = 'Kunz'

option.

As part of this new scheme we now ship a set of rates from NACREII `Xu et al 2013 <https://ui.adsabs.harvard.edu/abs/2013cgrs.conf..617X/abstract>`_. These rates do not, by default,
override the default NACRE rates. You must explicitly ask for them with ``filename_of_special_rate``.

There is now a new hook ``other_rate_get`` to provide a simple way to change an existing rate in a ``run_star_extras.f90``. Note this hook
only works on rates that are NOT currently in your rates_cache. It is recommended when using this option to set a custom ``rates_cache_dir`` otherwise the cache files in
``MESA_DIR`` will be over written.

The previous option:

::

    use_rate_3a = 'Fl87'

has been replaced with:

::

    use_3a_fl87 = .true.


net
---

There is a new hook ``other_net_derivs`` that allows for modifying the ``dydt`` term MESA computes for each zone inside ``net/``.
This allows adding changes in composition due to nuclear reactions that MESA could otherwise not handle or does not know about. 
This hook only works with soft networks (thus no ``approx`` nets). This hook requires many derivatives to be set, 
thus users should look at ``net_derivs.f90`` for reference to what needs setting.

There is now a hook ``other_split_burn`` for replacing MESA's split burn routine.

star
----

Diffusion coefficients for white dwarf interiors are now included based on
`Caplan et al. (2022) <https://ui.adsabs.harvard.edu/abs/2022MNRAS.tmpL..33C/abstract>`_.
By default, these coefficients are used for strong plasma coupling :math:`\Gamma > 10`, but there is an inlist option
to turn them off and revert to the previous default `Stanton & Murillo (2016) <https://ui.adsabs.harvard.edu/abs/2016PhRvE..93d3203S/abstract>`_
coefficients if desired.

Fixed a combination of bugs whereby the atmosphere data written to pulsation file formats (e.g. FGONG)
was incorrect or wouldn't work if ``tau_factor`` or ``atm_T_tau_opacity`` differed from their defaults
(``1.0`` and ``'fixed'``, respectively).

pgstar
------

Due to re-organization of the star_type derived type, all pgstar controls have been moved into a separate pgstar derived type.
If you access a pgstar option ``XX`` in your ``run_star_extras.f90`` then you need to replace ``s% XX`` with ``s% pg% XX``.

RSP
---

In r21.12.1 an experimental RSP solver feature was turned on by default, leading to convergence issues in nonlinear model integration. This is now turned off by default. Users that continue to use RSP in r.21.12.1 should include RSP_do_check_omega = .true. in the &controls section of their inlists to get rid of this issue.


Changes in r21.12.1
===================

.. _Backwards-incompatible changes r21.12.1:

Backwards-incompatible changes
------------------------------


.. note::

   A large amount of internal clean up has occurred since the last release.  This lists some of the most important changes, but the list is not exhaustive.



Simplification of energy equation options
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The desired form of the MESA energy equation is now selected via the control ``energy_eqn_option``.  The available options are
``'dedt'`` (default) and ``'eps_grav'``.  See the documentation at :ref:`reference/controls:energy_eqn_option` for more information about these forms.

The controls ``use_dedt_form_of_energy_eqn``, ``always_use_dedt_form_of_energy_eqn``, and ``use_eps_grav_form_of_energy_eqn`` were removed and replaced by the functionality of ``energy_eqn_option``.

Simplifications to the energy equation code mean that this selection applies globally (i.e., to all cells in the model and at all timesteps).

* The per-cell energy equation controls ``max_eta_for_dedt_form_of_energy_eqn`` and ``max_gamma_for_dedt_form_of_energy_eqn`` were removed.

* The form-switching control ``steps_before_always_use_dedt_form_of_energy_eqn`` was removed.


Name changes
~~~~~~~~~~~~

* The ``star_job`` option ``saved_model_name`` has been replaced with ``load_model_filename`` everywhere.

* The ``controls`` options ``power_c_burn_{lower,upper}_limit`` were replaced with the more generic ``power_z_burn_{lower,upper}_limit``.

* The ``controls`` option ``delta_lgL_phot_limit`` was renamed to ``delta_lgL_power_photo_limit`` ("phot" was easily confused with photosphere instead of photodisintegration).

* The ``controls`` options ``surf_w_div_w_crit_limit`` and ``surf_w_div_w_crit_tol`` were renamed to ``surf_omega_div_omega_crit_limit`` and ``surf_omega_div_omega_crit_tol``

* The core/layer mass values ``c_core_*``, ``c_rich_layer``, and
  ``o_core_*`` have been renamed to ``co_core_*``,
  ``co_rich_layer_*``, and ``one_core_*``.  This better reflects the
  typical carbon/oxygen and oxygen/neon compositions of these regions.
  This affects the names of both the relevant controls and history
  columns.

* The ``controls`` option ``use_d_eos_dxa`` was renamed to
  ``fix_d_eos_dxa_partials``.  This control originally had a broader
  function during the implementation of eos composition derivatives,
  but is now restricted to selecting whether we do a
  finite-difference-based fix up when on a component EOS that doesn't
  provide composition derivatives.

* The history and profile columns ``burn_*`` where replace with ``*_alpha``.

* History, profile, and binary history column files are now case insensitive. 


Removed options
~~~~~~~~~~~~~~~

* The time-smoothing scheme for mixing diffusion coefficients was removed.  All associated options (e.g., ``new_D_smooth_flag`` and ``D_smooth_replacement_fraction``) were removed.

* Removed option ``semiconvection_upper_limit_center_h1``. This can be implemented by setting ``s% alpha_semiconvection`` in ``run_star_extras.f90/extras_start_step``.

* Removed the option ``use_brunt_gradmuX_form``.  Alternative forms of the Brunt can be calculated using the ``other_brunt`` hook.

Removed history and profile columns
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A major clean up of the history and profile columns was undertaken.  Some of the removed values include:

* Removed profile columns ``total_energy`` and ``total_energy_integral``.


Relocation of eos hooks
~~~~~~~~~~~~~~~~~~~~~~~

The ``other_eos`` hooks have been removed from star.  See the ``eos`` section for information about their replacements.


Hook interface changes
~~~~~~~~~~~~~~~~~~~~~~

* The ``Teff`` argument has been removed from the ``other_surface_PT`` hook. (``Teff`` is instead available in the ``star_info`` pointer.)

* ``other_mesh_delta_coeff_factor`` no longer takes ``eps_h``, ``eps_he`` or ``eps_z`` as arguments.


Auto diff
~~~~~~~~~

We now make more extensive use of the new ``autodiff`` module for automatically differentiating variables. If you are using a hook
in your ``run_star_extras.f90`` then you will need to add ``use auto_diff`` to the top of your  ``run_star_extras.f90`` file.

If you see errors such as:

.. code-block:: fortran
  
  Error: Cannot convert REAL(8) to TYPE(auto_diff_real_star_order1) at (1)


Then this means you are missing the ``use auto_diff`` statement.

An example of using ``autodiff`` in a hook can be found at :ref:`autodiff example`


.. _Module-level changes r21.12.1:

Module-level changes
--------------------

astero
~~~~~~

Many of the one-dimensional arrays of mode data (e.g. ``l0_obs``) have
been consolidated into two-dimensional arrays (e.g. ``freq_target``)
in which the first index is the angular degree ``l``.  The following
controls in ``&astero_search_controls`` have changed:

+-----------------------+-----------------------+
+ Old                   + New                   +
+=======================+=======================+
+                       +                       +
+ ``nl0``               + ``nl(0)``             +
+                       +                       +
+ ``nl1``               + ``nl(1)``             +
+                       +                       +
+ ``nl2``               + ``nl(2)``             +
+                       +                       +
+ ``nl3``               + ``nl(3)``             +
+                       +                       +
+-----------------------+-----------------------+
+                       +                       +
+ ``l0_obs(:)``         + ``freq_target(0,:)``  +
+                       +                       +
+ ``l1_obs(:)``         + ``freq_target(1,:)``  +
+                       +                       +
+ ``l2_obs(:)``         + ``freq_target(2,:)``  +
+                       +                       +
+ ``l3_obs(:)``         + ``freq_target(3,:)``  +
+                       +                       +
+-----------------------+-----------------------+
+                       +                       +
+ ``l0_obs_sigma(:)``   + ``freq_sigma(0,:)``   +
+                       +                       +
+ ``l1_obs_sigma(:)``   + ``freq_sigma(1,:)``   +
+                       +                       +
+ ``l2_obs_sigma(:)``   + ``freq_sigma(2,:)``   +
+                       +                       +
+ ``l3_obs_sigma(:)``   + ``freq_sigma(3,:)``   +
+                       +                       +
+-----------------------+-----------------------+
+                       +                       +
+ ``iscan_factor_l0``   + ``iscan_factor(0)``   +
+                       +                       +
+ ``iscan_factor_l1``   + ``iscan_factor(1)``   +
+                       +                       +
+ ``iscan_factor_l2``   + ``iscan_factor(2)``   +
+                       +                       +
+ ``iscan_factor_l3``   + ``iscan_factor(3)``   +
+                       +                       +
+-----------------------+-----------------------+

The call signatures to the surface correction subroutines have also
changed, generally from

::

    subroutine get_some_freq_corr(...,
          nl0, l0_obs, l0_sigma, l0_freq, l0_freq_corr, l0_inertia,
          nl1, l1_obs, l1_sigma, l1_freq, l1_freq_corr, l1_inertia,
          nl2, l2_obs, l2_sigma, l2_freq, l2_freq_corr, l2_inertia,
          nl3, l3_obs, l3_sigma, l3_freq, l3_freq_corr, l3_inertia)

to

::

    subroutine get_some_freq_corr(...,
          nl, obs, sigma, freq, freq_corr, inertia)


binary
~~~~~~

There are new hooks ``other_binary_photo_read`` and
``other_binary_photo_write``.  These allow the user to save/restore
values in ``run_binary_extras``.


eos
~~~

There are new module-level eos hooks (see ``eos/other``) that replace
the star-level eos hooks (previously in ``star/other``).  Usage of
these hooks is similar to hooks in star.  However, the relevant
procedure pointer is part of the ``EOS_General_Info`` structure and
not the ``star_info`` structure.  Therefore, in ``extras_controls``,
the procedure pointer statement should look like ``s% eos_rq %
other_eos_results => my_other_eos_results``.  The boolean option
``use_other_eos_results`` controlling whether to use the hook is part
of the ``eos`` namelist rather than ``controls``.  For the first
required argument ``handle``, pass ``s% eos_handle``.  This ensures
that the routine uses the same configuration options as other calls
from star to the eos module.

The hook ``other_eos_component`` allows the user to replace all or
part of the MESA EOS by providing a new component EOS and to control
the location of the blends between this and the other component EOSes.
It is controlled by the option ``use_other_eos_component``.  The
user-provided routine must return a complete set of EOS results.  This
EOS component has the highest priority in the blend.  This hook
should be used along with the hook ``other_eos_frac``, which defines
the region over to use ``other_eos_component``.

The hook ``other_eos_results`` allows the user to modify the results
returned by the EOS.  The user-provided routine receives the results
from the EOS right before they are returned, after all components have
been evaluated.  This allows the user make minor modifications to the
results from the existing EOS without having to provide a full replacement.

Two alternative eos module entry points (``eosDT_HELMEOS_get`` and
``eosDT_ideal_gas_get``) and the star options that replaced the
standard eosDT calls to be with these routines
(``use_eosDT_ideal_gas`` and ``use_eosDT_HELMEOS``).  This enables
significant simplifications of eos_support.  Restriction to a single
component EOS can be achieved through the eos namelist options and
replacement of the EOS should be performed through the other hook.

The HELM table was updated to a new, larger 100 points per decade
version.

The HELM-related controls ``logT_ion_HELM``, ``logT_neutral_HELM``, and
``max_logRho_neutral_HELM`` were removed.  These were used in an
now-unsupported variant of HELM that blended the normal, fully-ionized
HELM and a neutral version (which dropped the electron-positron terms).

The HELM-related controls ``always_skip_elec_pos`` and
``always_include_elec_pos`` were combined in the
simplified control ``include_elec_pos`` which defaults to ``.true.``.

There is a new backstop EOS (``ideal``) which analytically models an ideal ion gas with radiation pressure.
The purpose of this EOS is to provide coverage over the whole density-temperature plane for times when MESA needs
to run to extreme densities or temperatures.
No electrons are included in this EOS.


kap
~~~

The call signatures of ``kap_get`` and the hook ``other_kap_get`` have
changed.  The set of arguments is now conceptually equivalent between
the two subroutines.  The inputs include the density, temperature, and
full composition vector.  The free electron/positron number and the
electron degeneracy parameter (and their derivatives) are also
required.  The outputs include the opacity and its derivatives as well
as information about the fractions of various opacity sources used in
the blended opacity.

::

      subroutine kap_get( &
         handle, species, chem_id, net_iso, xa, &
         logRho, logT, &
         lnfree_e, d_lnfree_e_dlnRho, d_lnfree_e_dlnT, &
         eta, d_eta_dlnRho, d_eta_dlnT , &
         kap_fracs, kap, dlnkap_dlnRho, dlnkap_dlnT, dlnkap_dxa, ierr)

         ! INPUT
         integer, intent(in) :: handle ! from alloc_kap_handle; in star, pass s% kap_handle
         integer, intent(in) :: species
         integer, pointer :: chem_id(:) ! maps species to chem id
         integer, pointer :: net_iso(:) ! maps chem id to species number
         real(dp), intent(in) :: xa(:) ! mass fractions
         real(dp), intent(in) :: logRho ! density
         real(dp), intent(in) :: logT ! temperature
         real(dp), intent(in) :: lnfree_e, d_lnfree_e_dlnRho, d_lnfree_e_dlnT
            ! free_e := total combined number per nucleon of free electrons and positrons
         real(dp), intent(in)  :: eta, d_eta_dlnRho, d_eta_dlnT
            ! eta := electron degeneracy parameter from eos

         ! OUTPUT
         real(dp), intent(out) :: kap_fracs(num_kap_fracs)
         real(dp), intent(out) :: kap ! opacity
         real(dp), intent(out) :: dlnkap_dlnRho ! partial derivative at constant T
         real(dp), intent(out) :: dlnkap_dlnT   ! partial derivative at constant Rho
         real(dp), intent(out) :: dlnkap_dxa(:) ! partial derivative w.r.t. species
         integer, intent(out) :: ierr ! 0 means AOK.


The Compton scattering opacity routine has been updated to use the prescription of
`Poutanen (2017) <https://ui.adsabs.harvard.edu/abs/2017ApJ...835..119P/abstract>`_.

The conductive opacity routine has been updated to include the corrections from 
`Blouin et al. (2020) <https://ui.adsabs.harvard.edu/abs/2020ApJ...899...46B/abstract>`_
for H and He in the regime of moderate coupling and moderate degeneracy.
These are on by default, controlled by the kap option ``use_blouin_conductive_opacities``.


There are new module-level kap hooks (see ``kap/other``) that allow
individual components of the opacity module to be replaced with a
user-specified routine given in run_star_extras.  Usage of these hooks
is similar to hooks in star.  However, the relevant procedure pointer
is part of the ``Kap_General_Info`` structure and not the
``star_info`` structure.  Therefore, in ``extras_controls``, the
procedure pointer statement should look like ``s% kap_rq %
other_elect_cond_opacity => my_routine``.  The boolean option
``use_other_elect_cond_opacity`` controlling whether to use the hook
is part of the ``kap`` namelist rather than ``controls``.  For the
first required argument ``handle``, pass ``s% kap_handle``.  This
ensures that the routine uses the same configuration options as other
calls from star to the kap module.


neu
~~~

The call signature of other_neu has changed. You no longer need to pass in z2bar


The value of the Weinberg angle was updated to be be consistent with CODATA 2018.


net
~~~

The screening mode ``classic_screening`` has been removed. Anyone using other_net_get needs
to remove ``theta_e_for_graboske_et_al`` from its argument list.

The options ``reuse_rate_raw`` and  ``reuse_rate_screened`` have been removed from other_net_get (and eval_net)


rates
~~~~~

The format for custom weak rate tables (see e.g., ``data/rates_data/rate_tables/weak_rate_list.txt`` and test suite case ``custom_rates``) no longer supports the (previously optional) Coulomb correction datasets ``delta_Q`` and ``Vs``.

When this capability was first added, the energetics associated with
the change in the composition were calculated in ``rates`` and
included in ``eps_nuc``.  This meant the ``rates`` module needed to
have access to information about the Coulomb-induced shifts in the
electron and ion chemical potentials.

After the changes in the definition of ``eps_nuc`` and the energy
equation described in |MESA V|, the energetics associated with the
changing composition are self-consistently accounted for in the energy
equation using information provided by the MESA EOS.  Therefore, the
ability to provide these unneeded and unused quantities has been
removed.


.. _Other changes r21.12.1:

Other changes
-------------

* Analogous to ``kap_frac_Type2``, information about the fractional
  contribution of the lowT tables, highT tables, and Compton opacities
  to the final result from the opacity module are now included in
  star_info arrays and profile columns with the names
  ``kap_frac_lowT``, ``kap_frac_highT``, ``kap_frac_Compton``.

* The control ``format_for_FGONG_data`` has been replaced by the
  integer ``fgong_ivers``, which can be either 300 or 1300.  This
  enforces adherence to the FGONG standard.  In addition, users can
  now set the four-line plain-text header of FGONG output using the
  new controls ``fgong_header(1:4)``.

* ``mixing_type`` now reports the mixing process that generates the
  largest D_mix, rather than prioritizing convection and thermohaline
  mixing over all others.

* Added profile panel and history panel controls in pgstar to specify
  same yaxis range for both left and right axes (e.g.,
  Profile_Panels1_same_yaxis_range(1) = .true.)

* Experimental options have been moved into ``*_dev.defaults`` files
  and experimental test cases are now prefixed with ``dev_``.  These
  options and test cases are not ready for general use.

* The ``ionization`` module has been removed.
  The ``eval_typical_charge`` routine has been moved into
  ``mod_typical_charge.f90`` within the ``star`` module.
  The ``eval_ionization`` routine is no longer
  supported, as it was untested, undocumented, and unused.

* A new module ``hdf5io`` for working with HDF5 files has been added.

* The controls ``diffusion_gamma_full_{on,off}`` are no longer used by
  default.  The EOS now returns phase information and by default that
  EOS phase will automatically turn off diffusion for crystallized
  material.

* The `issue with the value of free_e when using FreeEOS <https://lists.mesastar.org/pipermail/mesa-users/2021-February/012394.html>`__ has been corrected.  Thanks to Jason Wright for the report.

* An ``other_screening`` hook was added.

* All parts of test suite cases are now run by default.  To skip
  running the optional inlists, set the environment variable
  ``MESA_SKIP_OPTIONAL`` (to any value).  Previously, optional parts
  were skipped by default, and running all parts required setting
  ``MESA_RUN_OPTIONAL``.

* The headers for history and profile data now contain the value of Msun (grams), Rsun (cm), and Lsun (erg/s) used.

* A bug has been identified and fixed in the ``Brown_Garaud_Stellmach_13``
  thermohaline mixing routine. The routine was meant to use
  Newton-Raphson relaxation to converge to a solution for the Nusselt
  number based on an initial guess from the asymptotic analysis in
  Appendix B of
  `Brown, Garaud, & Stellmach (2013) <https://ui.adsabs.harvard.edu/abs/2013ApJ...768...34B>`_.
  However, a bug previously caused the routine to immediately return the
  asymptotic guess and skip the NR relaxation step. The asymptotic
  guess is usually fairly accurate, so this usually still produced a
  thermohaline result that was fairly close to the right answer, but
  the bug has been fixed now so that the NR relaxation is applied as
  well.


Changes in r15140
=================

.. _Backwards-incompatible changes r15140:

Backwards-incompatible changes
------------------------------

Addition of eos and kap namelists
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The options associated with the ``eos`` and ``kap`` modules have been
moved into their own namelists.  (That is, there now exist ``&eos``
and ``&kap`` at the same level as ``&star_job`` and ``&controls``.)
User inlists will need to be updated.  See :ref:`Module-level changes r15140`
for more specific information.

If you previously accessed the values of eos/kap related options from
``star_job`` or ``controls`` via run_star_extras, you will need to
adjust your code to access the option values using the pointers to the
``EoS_General_Info`` and ``Kap_General_Info`` structures.  These are
exposed in star as ``s% eos_rq`` and ``s% kap_rq``, respectively.  So
for example, the inlist value of ``Zbase`` is now accessible via ``s% kap_rq% Zbase``
(instead of ``s% Zbase``).


Some file suffixes changed to .f90
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Many source file names have been changed to have an .f90 suffix.  For
users, the most important changes are to the star and binary work
directories.

In an existing star work directory (i.e., a copy of ``star/work`` or
star test suite case), rename the files

+ ``src/run.f`` → ``src/run.f90``
+ ``src/run_star_extras.f`` → ``src/run_star_extras.f90``

In an existing binary work directory (i.e., a copy of
``binary/work`` or binary test suite case), rename the files

+ ``src/binary_run.f`` → ``src/binary_run.f90``
+ ``src/run_star_extras.f`` → ``src/run_star_extras.f90``
+ ``src/run_binary_extras.f`` → ``src/run_binary_extras.f90``

Changes to local makefiles that are not part of MESA might also need
to be updated to reflect these changes.


Removal of backups
~~~~~~~~~~~~~~~~~~

MESA no longer has the concept of a "backup".  (In a backup, after the
failure of a retry, MESA would return to the previous model and evolve
it with a smaller timestep.)

Models that previously relied on the use of backups in order to
complete should instead use appropriate timestep controls such that
retries alone are sufficient to enable the model to run.

All backup-related options and output quantities have been removed.
Users migrating inlists or ``history_column.list`` files from previous
MESA versions will need to remove these options, all of which contain
the string "backup".


Changes to solver reporting
~~~~~~~~~~~~~~~~~~~~~~~~~~~

MESA can report information about the progress of the iterative
Newton–Raphson solution process that forms a key part of taking a
timestep.  The names of numerous options related to the solver have
changed.  These changes follow two main patterns.

First, the word "newton" was replaced with the word "solver".  For
example, the history column that records the number of iterations
changed from ``num_newton_iterations`` to ``num_solver_iterations``.
The controls option that defines a number iterations above which to
reduce the timestep changed from ``newton_iterations_limit`` to
``solver_iters_timestep_limit`` and the terminal output correspondingly
shows the message ``solver iters`` instead of ``newton iters``.  (The
control ``newton_iterations_hard_limit`` was removed and not renamed.)

Second, the word "hydro" was removed or replaced with the word
"solver" in the controls related to monitoring the solver internals.
For example, the control ``report_hydro_solver_progress`` is now
``report_solver_progress`` and ``report_hydro_dt_info`` is now
``report_solver_dt_info``.  The use of these and other related
controls is described :ref:`in the developer documentation <developing/debugging:Diagnosing Solver Struggles>`.



Changes to eps_grav and eps_mdot
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A new method for handling the energetics associated with mass changes
in MESA models was presented in |MESA V|, Section 3.2.  The approach
discussed therein, incorporated in a term named ``eps_mdot``, has now
become standard.  As such, the option ``use_eps_mdot`` has been
removed (because it is now effectively always true).

This ``eps_mdot`` approach supersedes the approach described in
|MESA III|, Section 7, and so that implementation has been removed.  This
resulted in the removal of the ``&controls`` options

+ ``eps_grav_time_deriv_separation``
+ ``zero_eps_grav_in_just_added_material``
+ ``min_dxm_Eulerian_div_dxm_removed``
+ ``min_dxm_Eulerian_div_dxm_added``
+ ``min_cells_for_Eulerian_to_Lagrangian_transition``
+ ``fix_eps_grav_transition_to_grid``

the history columns
  
+ ``k_below_Eulerian_eps_grav``
+ ``q_below_Eulerian_eps_grav``
+ ``logxq_below_Eulerian_eps_grav``
+ ``k_Lagrangian_eps_grav``
+ ``q_Lagrangian_eps_grav``
+ ``logxq_Lagrangian_eps_grav``

and the profile columns
  
+ ``eps_grav_h_effective``
+ ``eps_mdot_sub_eps_grav_h_effective``
+ ``eps_mdot_rel_diff_eps_grav_h_effective``
+ ``eps_grav_h``
+ ``eps_mdot_sub_eps_grav_h``
+ ``eps_mdot_rel_diff_eps_grav_h``


Removal of lnPgas_flag
~~~~~~~~~~~~~~~~~~~~~~

The option to use gas pressure instead of density as a structure
variable has been removed.  Users migrating inlists from previous MESA
versions will need to remove these options, all of which contain the
string "lnPgas_flag".


Removal of logQ limits
~~~~~~~~~~~~~~~~~~~~~~

As a consequence of the changes to ``eos``, ``star`` no longer
enforces limits on the quantity logQ (``logQ = logRho - 2*logT + 12`` in cgs).
Therefore the ``controls`` options

- ``logQ_limit``
- ``logQ_min_limit``

and the ``pgstar`` option

- ``show_TRho_Profile_logQ_limit``

have been removed.

The removal of these controls does not indicate that the EOS is
reliable at all values of logQ.  Users should consult :ref:`the
description of the component EOSes and the regions in which they are
applied <eos/overview:Overview of eos module>` to understand if MESA provides
a suitable EOS for the conditions of interest.


Removal of GR factors
~~~~~~~~~~~~~~~~~~~~~

The control ``use_gr_factors`` and corresponding code has been
removed.  (This provided only a simple correction to the momentum
equation and not a full GR treatment of the stellar structure
equations.)  Users wishing to include GR corrections to MESA's
Newtonian equations can achieve the same effect by using the
``other_cgrav`` or ``other_momentum`` hooks.  For an example, see the
neutron star test cases (``ns_h``, ``ns_he``, and ``ns_c``).


Change in STELLA file output
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The options to create output files suitable for input to STELLA have
been removed from ``MESA/star`` and the ``star_job`` namelist.  These
capabilities are now included as part of the ``ccsn_IIp`` test case
(see ``inlist_stella`` and ``run_star_extras.f90``).  Users desiring
STELLA-format output should re-use the code from that example.

This affects the options

- ``save_stella_data_for_model_number``
- ``save_stella_data_when_terminate``
- ``save_stella_data_filename``
- ``stella_num_points``
- ``stella_nz_extra``
- ``stella_min_surf_logRho``
- ``stella_min_velocity``
- ``stella_skip_inner_dm``
- ``stella_skip_inner_v_limit``
- ``stella_mdot_years_for_wind``
- ``stella_mdot_for_wind``
- ``stella_v_wind``
- ``stella_show_headers``

  
Removal of mesh adjustment parameters around convective boundaries
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Controls matching the following patterns, which adjust the mesh
resolution around convective boundaries, have been removed:

- ``xtra_coef_czb_full_{on,off}``
- ``xtra_coef_{a,b}_{l,u}_{n,h,he,z}b_czb``
- ``xtra_dist_{a,b}_{l,u}_{n,h,he,z}b_czb``
- ``xtra_coef_scz_above_{n,h,he,z}b_cz``

Convective boundaries can be resolved using a custom mesh-spacing
function or ``mesh_delta`` coefficient.  The
``simplex_solar_calibration`` test case has an example custom
mesh-spacing function.


Change to ``mixing_type`` codes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ``mixing_type`` codes (defined in ``const/public/const_def.f90``)
have changed.  User code and/or analysis routines (e.g., scripts
interpreting the ``mixing_type`` profile column) may need to be
revised.  We recommend that users use the ``mixing_type`` variables
rather than the corresponding integers in their own code. e.g. rather
than writing
::

    if (mixing_type == 1) then

use
::

    if (mixing_type == convective_mixing) then

assuming ``use const_def`` appears somewhere, as in the default
``run_star_extras.f90``.

Limitations on use of varcontrol_target
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A new variable ``min_allowed_varcontrol_target`` (default 1d-4) has
been introduced to discourage the use of small values of
``varcontrol_target``.  MESA will exit with an error if the value is
below this threshold.

The value of ``varcontrol`` is an unweighted average over all cells of the
relative changes in the structure variables.  For situations that need
tighter timestep limits, there are many specific timestep controls
that should be used instead of reducing the general target.  The use
of controls that specifically apply to the problem being studied will
typically provide more effective and efficient timestep limiters.  In
addition, small values of ``varcontrol_target`` can lead to poor
performance when it forces the size of the step-to-step corrections to
become too small.

The option ``varcontrol_target`` is NOT the recommended way to push
time resolution to convergence levels. To perform temporal convergence
studies, use the new control ``time_delta_coeff``, which acts as a
multiplier for timestep limits (analogous to ``mesh_delta_coeff`` for
spatial resolution).

One strategy for choosing effective timestep limits is to first set
``varcontrol_target = 1d-3``.  Then add some additional specific
timestep limits relevant to the problem.  Do a run, watching the
reason for the timestep limits and the number of retries.  Summary
information about the conditions that limited the timestep can be
printed at the end of run using the ``star_job`` option
``show_timestep_limit_counts_when_terminate``.  Repeat the runs,
adding/removing or adjusting timestep limits until there are few
retries and few places where the timestep is limited by varcontrol.
Finally, repeat the calculation with a smaller value of
``time_delta_coeff`` (e.g., 0.5) and compare the results to gain
confidence that they are numerically converged.

.. _Module-level changes r15140:

Module-level changes
--------------------

astero
~~~~~~

Material previously present in ``star/astero`` and test cases using
these capabilities have been promoted into their own module.

The ``csound_rms`` observational constraint has been removed.

The options for executing an arbitrary shell script
(``shell_script_num_string_char`` and
``shell_script_for_each_sample``) have been removed.  The usual use
for these options—renaming output files at the end of each sample—can
be replicated using the system tools available through
``utils_lib``.  For example, the following ``extras_after_evolve``
in ``run_star_extras.f90`` moves the best profile and FGONG file
to ``outputs/sample_#.{profile,fgong}``.
::

      subroutine extras_after_evolve(id, ierr)
         use astero_def
         use utils_lib, only: mv
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         character (len=256) :: format_string, num_string, basename
         ierr = 0

         write(format_string,'( "(i",i2.2,".",i2.2,")" )') num_digits, num_digits
         write(num_string,format_string) sample_number+1 ! sample number hasn't been incremented yet
         basename = trim(sample_results_prefix) // trim(num_string)
         call mv(best_model_fgong_filename, trim(basename) // trim('.fgong'), skip_errors=.true.)
         call mv(best_model_profile_filename, trim(basename) // trim('.profile'), skip_errors=.true.)
         
      end subroutine extras_after_evolve

turb
~~~~

This new module implements local theories of turbulence, including
MLT, TDC, semiconvection, and thermohaline turbulence. These used to be
a part of ``star``. TDC (which stands for time-dependent convection) is
now the recommended method for situations where the time dependence of
convection must be taken into account. Other methods for time dependent
convection present in the code have been removed, including the options
min_T_for_acceleration_limited_conv_velocity and set_conv_vel_flag.
TDC can be turned on with the option ``MLT_option = "TDC"`` in the 
``controls`` section of an inlist.

Users will not generally
need to interact with this module, but it can be used within
run_star_extras by writing ``use turb``.

auto_diff
~~~~~~~~~

This new module provides Fortran types that support algorithmic
differentiation via operator overloading. Users will not generally
need to interact with this module, but it can be used within
run_star_extras to make derivatives easier to calculate (e.g. in the
implicit hooks like ``other_surface``).

Usage is by writing ``use auto_diff``. This imports types such as
``auto_diff_real_4var_order1``, which supports first-order derivatives
with respect to up to four independent variables.
A variable of this type could be declared via::

    type(auto_diff_real_4var_order1) :: x

This variable then holds five fields: ``x%val`` stores the value of ``x``.
``x%d1val1`` stores the derivative of `x` with respect to the first independent
variable. ``x%d1val2`` is the same for the second independent variable, and so on.
All ``d1val_`` fields are initialized to zero when the variable is first set.

Once an auto_diff variable it initialized, all mathematical operations can be performed
as they would be on a ``real(dp)`` variable. auto_diff variables also interoperate with
``real(dp)`` and ``integer`` types.

So for instance in the following ``f%d1val1`` stores df/dx and ``f%d1val2`` stores df/dy.
::
   
    x = 3d0
    x%d1val1 = 1d0
    
    y = 2d0
    y%d1val2 = 1d0
    
    f = exp(x) * y + x + 4

Similar types are included supporting higher-order and mixed-partial
derivatives.  These derivatives are accessed via e.g. ``d2val1``
(d²f/dx²), ``d1val1_d2val2`` (d³f/dx dy²).


const
~~~~~

The ``const`` module has been updated to account for the revision of
the SI and now uses CODATA 2018 values of the physical constants.

For astronomical constants, MESA follows IAU recommendations. MESA
adopts nominal solar and planetary quantities from IAU 2015 Resolution
B3 and now follows the recommended procedure of deriving nominal solar
and planetary masses from the mass parameters :math:`(GM)` and the
adopted value of :math:`G`.

As a result of these changes, most constants now have slightly
different values than in previous MESA versions. For example, |Msun|
changed from 1.9892e33 g to 1.9884e33 g.


eos
~~~

EOS-related options have been moved into their own ``eos`` namelist.
The :ref:`module controls <eos/defaults:eos module controls>` and their default
values are contained in the file ``eos/defaults/eos.defaults``.

The PTEH EOS has been removed.  Tables from the FreeEOS project now
provide coverage of similar conditions.

The region covered by the PC EOS has been increased.  The boundaries
of the region where PC is used no longer consider composition and so
now include H/He-dominated material.  The upper limit of the region
where PC is used is now determined using the electron Coulomb coupling
parameter and generally corresponds to higher temperatures than the
previous approach.

For more information about the component EOSes and the regions in
which they are applied, see the :ref:`new overview of the EOS module
<eos/overview:Overview of eos module>`.


gyre
~~~~

GYRE has been upgraded to version 6.0.  See the `GYRE Documentation
<https://gyre.readthedocs.io/en/latest/index.html>`__ for information
about this release.

kap
~~~

Opacity-related options have been moved into their own ``kap`` namelist.
The :ref:`module controls <kap/defaults:kap module controls>` and their default
values are contained in the file ``kap/defaults/kap.defaults``.


The OPAL Type 2 opacity tables are now on by default
(``use_Type2_opacities = .true.``).  These tables separately account
for carbon and oxygen enhancements.  Since this is especially
important during core helium burning, the default transition from the
OPAL Type 1 tables to the Type 2 tables occurs when material becomes
nearly hydrogen free.  As a result of this change, by default, users
are required to specify the base metallicity of material using the
``kap`` namelist control ``Zbase``.  Usually, this physically
corresponds to the initial metallicity of the star.


For more information about the opacity tables and how they are
combined, see the :ref:`new overview of the kap module <kap/overview:Overview of
kap module>`.

rates & net
~~~~~~~~~~~

A number of rates have had their defaults switched to using JINA's REACLIB.

When using a custom user rate (i.e from a rate table) the reverse rate is now computed in detailed
balance from the user rate. Previously the reverse rate was computed using the default rate choice.

A bug with burning li7 at low temperatures rate has been fixed. Users stuck using previous versions of MESA and 
a soft network (something that is not an approx net) should add these lines to their nuclear network as a fix until they
can update to a newer MESA:
::

    remove_reaction(r_h1_li7_to_he4_he4)
    add_reaction(r_li7_pa_he4)

With thanks to Ian Foley for the bug report.

We now define the forward reaction to always be the exothermic reaction, not the reaction as defined by REACLIB.
This fixes an issue with exothermic photo-disintegrations which would generate wrong values when computed
in detailed balance.

A lot of work has been done getting operator split burning (op_split_burn = .true.) to work.
This option can provide a large speed up during advanced nuclear burning stages. See the various split_burn
test cases for examples.

.. _Other changes r15140:

Other changes
-------------

* Saved model files now contain a ``version_number`` field in their
  header.  This indicates the version of MESA that was used to
  generate the model.

* binary now automatically writes photo (restart) files at the end of
  the run.

* If not provided with an argument, the binary ``./re`` script will
  now restart from the most recent photo (determined by filesystem
  modification time).  The star ``./re`` script also has this behavior
  as of r12778.

* The test case for building C/O white dwarf models has been
  overhauled to be more robust. See documentation for the new version
  in :ref:`make_co_wd`.

* The builder for NS envelopes (test case ``neutron_star_envelope``)
  has been replaced with a more general envelope builder (test case
  ``make_env``).  The test cases ``ns_{h,he,c}`` have been overhauled
  to start from these new models.

* Added ``other_remove_surface``. This routine is called at the start
  of a step and returns an integer k. All cells with j < k will be removed
  from the model at the start of the step, making cell k the new surface.

* Installations are now blocked from using sudo. This is generally not what you want to
  use to fix installation issues. If you want to install MESA in a root location
  then you will need to edit out the check in the install file.

* The install script now blocks attempts to use a ``MESA_DIR`` which contains spaces in it.
  This has never really worked as makefiles can not handle the spaces. To work round this
  either move ``MESA_DIR`` to a folder location with no spaces in its path or symlink
  your ``MESA_DIR`` to another location with no spaces in its path and set ``MESA_DIR``
  to point at the symlink.

* The option to create a pre main sequence model now relaxes the model until
  a radiative  core forms. This is activated with the ``star_job`` option
  ``pre_ms_relax_to_start_radiative_core``, which can be set to .false. to
  restore the old behavior.

.. _Acknowledgments r15140:

Acknowledgments
---------------

Thanks to all who reported problems and asked or answered questions on
mesa-users.  Special thanks to Siemen Burssens, Mathias Michielsen,
Joey Mombarg, Mathieu Renzo, and Samantha Wu for their assistance in
testing pre-release versions.


Changes in r12778
=================

This section describes changes that occurred since r12115.

SDK changes (Version 20.3.1 or later required)
----------------------------------------------

To use the this MESA release, you must upgrade your SDK to 20.3.1.

In previous releases of MESA, we have included the `CR-LIBM library <https://hal-ens-lyon.archives-ouvertes.fr/ensl-01529804/file/crlibm.pdf>`__
to provide versions of standard math functions (exp, log, sin, etc)
that guarantee correct rounding of floating-point numbers. In this new
release, we made the decision to move CR-LIBM into the software
development kit (SDK), where it properly belongs and can be maintained
as one of the pre-requisites of MESA.

This means that to use this release (and subsequent releases) of MESA,
you'll need to upgrade to version 20.3.1 of the SDK or later. MESA
checks the SDK version during compilation, and will stop with an error
message if the SDK is too old.

.. _Backwards-incompatible changes r12278:

Backwards-incompatible changes
------------------------------

Replacement of crlibm_lib with math_lib
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

MESA now talks to CR-LIBM via an intermediate module called
``math_lib``. To make sure any code you add can properly access the
CR-LIBM math routines, you'll need to make sure that a ``use
math_lib`` statement appears in the preamble of the file. At the same
time, you should remove any ``use crlibm_lib`` statements, as they will no
longer work (and are not needed).  With ``math_lib``, the names of the
correctly rounded math functions are the same as the Fortran
intrinsics (i.e., they no longer have a ``_cr`` suffix).

Existing ``run_star_extras``, ``run_binary_extras``, or other
user-written code will need to be updated.  To migrate, you should
replace ``use crlibm_lib`` with ``use math_lib`` and remove the ``_cr``
suffix from any math functions (e.g., ``exp_cr`` → ``exp``).


Removal of DT2 and ELM EOS options
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ELM and DT2 EOS options have been removed.  (These options were
underpinned by HELM and OPAL/SCVH data, but used bicubic spline
interpolation in tables of lnPgas, lnS, and lnE as a way to get
numerically accurate 1st and 2nd partial derivatives with respect to
lnRho and lnT.  A more detailed description can be found in previous
release notes and Appendix A.1 of |MESA V|.) These options were
introduced in r10398 and were turned on by default in r11532.

The numerical issues that ELM/DT2 were designed to address have been
dealt with via another approach.  MESA now separately treats quantities
that appear in the equations (and happen to be partials) and the
places where these theoretically equivalent, but numerically different
quantities appear in the Jacobian (as partials of other quantities
that appear in the equations).  This is an implementation detail that
should be transparent to users.

This change has two pleasant side effects.  One, it lowers the memory
demands of many MESA models, which should aid users of virtualized,
containerized, or otherwise memory-constrained systems.  Two, it
removes small, interpolation-related wiggles that were present in
partial derivative quantities such as :math:`\Gamma_1`.

These changes may require inlists that made use of DT2/ELM related
options to be updated.

The following ``controls`` options have been deleted:

  * ``use_eosDT2``
  * ``max_logT_for_max_logQ_eosDT2``
  * ``max_logQ_for_use_eosDT2``

  * ``use_eosELM``
  * ``logT_max_for_ELM``
  * ``logQ_min_for_ELM``
  * ``check_elm_abar_zbar``
  * ``check_elm_helm_agreement``


The following ``star_job`` options have been renamed:

  * ``eosDT2PTEH_use_linear_interp_for_X`` to ``eosPTEH_use_linear_interp_for_X``
  
The following ``controls`` options have been renamed/removed, as well
as moved to ``star_job`` (see next entry):

  * ``logRho_max_for_all_PTEH_or_DT2`` to ``logRho_max_for_all_PTEH``
  * ``logRho_max_for_any_PTEH_or_DT2`` to ``logRho_max_for_any_PTEH``
  * ``logQ_max_for_low_Z_PTEH_or_DT2`` (removed)
  * ``logQ_max_for_high_Z_PTEH_or_DT2`` to ``logQ_max_for_PTEH``


Change in location of PTEH EOS options
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Options that modify the parameters associated with the PTEH EOS have
be moved from ``controls`` to ``star_job``.  This brings PTEH in line
with the behavior of the other component EOSes.

If you explicitly set any of following options in your inlist, you
will need to move them from ``controls`` to ``star_job``.  Their
meaning and default values remain unchanged.

  * ``use_eosPTEH_for_low_density``
  * ``use_eosPTEH_for_high_Z``
  * ``Z_for_all_PTEH``
  * ``Z_for_any_PTEH``
  * ``logRho_min_for_all_OPAL``
  * ``logRho_min_for_any_OPAL``
  * ``logRho_max_for_all_PTEH``
  * ``logRho_max_for_any_PTEH``

In addition, you must add the new option ``set_eosPTEH_parameters =
.true.`` to ``star_job`` to indicate that these values should override
the eos module-level defaults.

The removal of DT2 (see previous entry) has also resulted in the
change that the ``controls`` option ``logQ_max_for_low_Z_PTEH_or_DT2`` has
been removed (as it applied primarily to DT2) and
``logQ_max_for_high_Z_PTEH_or_DT2`` (which applied primarily to PTEH)
has been renamed to ``logQ_max_for_PTEH`` and moved from ``controls``
to ``star_job``.


New overshooting controls
~~~~~~~~~~~~~~~~~~~~~~~~~

The new controls for overshooting, briefly described in the release notes of version 12115, are now the default in MESA (and the old controls have been removed). All test_suite cases now use these new controls.

There are two schemes implemented in MESA to treat overshooting: a step overshoot scheme and an exponential scheme that follows
`Herwig (2000) <https://ui.adsabs.harvard.edu/abs/2000A%26A...360..952H/abstract>`__.

The old "double exponential overshoot scheme" is no longer accessible through simple controls.  An example of how to implement such a scheme via the ``other_overshooting_scheme`` hook is contained in the ``other_physics_hooks`` test suite case.

The new overshooting controls are based on convection-zone and convection-boundary matching criteria.
In the new set of controls, for each convective boundary it is possible
to define an ``overshoot_zone_type``, ``overshoot_zone_loc`` and an
``overshoot_bdy_loc``, as well as values for the overshooting parameters.

The permitted values are the following:

  * ``overshoot_scheme = exponential, step``
  * ``overshoot_zone_type = burn_H, burn_He, burn_Z, nonburn, any``
  * ``overshoot_zone_loc = core, shell, any``
  * ``overshoot_bdy_loc = bottom, top, any``

The following controls assign values for the diffusive or step
overshooting parameters:

  * ``overshoot_f``
  * ``overshoot_f0``
  * ``overshoot_D0``
  * ``overshoot_Delta0``

overshoot_f0 is defined so that the switch from convective mixing to overshooting happens at a distance overshoot_f0*Hp into the convection zone from the estimated location where `grad_ad = grad_rad`, where Hp is the pressure scale height at that location.

For exponential overshoot, D(dr) = D0*exp(-2*dr/(overshoot_f*Hp0) where D0 is the diffusion coefficient D at point r0, Hp0 is the scale height at r0.

For step overshoot:
overshooting extends a distance overshoot_f*Hp0 from r0 with constant diffusion coefficient  D = overshoot_D0 + overshoot_Delta0*D_ob
where D_ob is the convective diffusivity at the bottom (top) of the step overshoot region for outward (inward) overshooting.

These "new" controls replace the following "old" controls:

  * ``overshoot_f_above_nonburn_core``
  * ``overshoot_f0_above_nonburn_core``
  * ``overshoot_f_above_nonburn_shell``
  * ``overshoot_f0_above_nonburn_shell``
  * ``overshoot_f_below_nonburn_shell``
  * ``overshoot_f0_below_nonburn_shell``
  * ``overshoot_f_above_burn_h_core``
  * ``overshoot_f0_above_burn_h_core``
  * ``overshoot_f_above_burn_h_shell``
  * ``overshoot_f0_above_burn_h_shell``
  * ``overshoot_f_below_burn_h_shell``
  * ``overshoot_f0_below_burn_h_shell``
  * ``overshoot_f_above_burn_he_core``
  * ``overshoot_f0_above_burn_he_core``
  * ``overshoot_f_above_burn_he_shell``
  * ``overshoot_f0_above_burn_he_shell``
  * ``overshoot_f_below_burn_he_shell``
  * ``overshoot_f0_below_burn_he_shell``
  * ``overshoot_f_above_burn_z_core``
  * ``overshoot_f0_above_burn_z_core``
  * ``overshoot_f_above_burn_z_shell``
  * ``overshoot_f0_above_burn_z_shell``
  * ``overshoot_f_below_burn_z_shell``
  * ``overshoot_f0_below_burn_z_shell``
  * ``step_overshoot_f_above_nonburn_core``
  * ``step_overshoot_f_above_nonburn_shell``
  * ``step_overshoot_f_below_nonburn_shell``
  * ``step_overshoot_f_above_burn_h_core``
  * ``step_overshoot_f_above_burn_h_shell``
  * ``step_overshoot_f_below_burn_h_shell``
  * ``step_overshoot_f_above_burn_he_core``
  * ``step_overshoot_f_above_burn_he_shell``
  * ``step_overshoot_f_below_burn_he_shell``
  * ``step_overshoot_f_above_burn_z_core``
  * ``step_overshoot_f_above_burn_z_shell``
  * ``step_overshoot_f_below_burn_z_shell``
  * ``step_overshoot_D``
  * ``step_overshoot_D0_coeff``

   
The "new" control ``overshoot_D_min`` replaces the "old"  control  ``D_mix_ov_limit``.

The "new" control ``overshoot_brunt_B_max`` replaces the "old"  control  ``max_brunt_B_for_overshoot``.

The "new" control ``overshoot_mass_full_on`` replaces the "old"  control  ``mass_for_overshoot_full_on``.

The "new" control ``overshoot_mass_full_off`` replaces the "old"  control  ``mass_for_overshoot_full_off``.

The following example will apply exponential overshoot, with f = 0.128
and f0 = 0.100, at the bottom of non-burning convective shells; and
exponential overshoot, with f = 0.014 and f0 = 0.004, at all other
convective boundaries.

::

  overshoot_scheme(1) = 'exponential'
  overshoot_zone_type(1) = 'nonburn'
  overshoot_zone_loc(1) = 'shell'
  overshoot_bdy_loc(1) = 'bottom'
  overshoot_f(1) = 0.128
  overshoot_f0(1) = 0.100

  overshoot_scheme(2) = 'exponential'
  overshoot_zone_type(2) = 'any'
  overshoot_zone_loc(2) = 'any'
  overshoot_bdy_loc(2) = 'any'
  overshoot_f(2) = 0.014
  overshoot_f0(2) = 0.004

Other examples are illustrated in the test_suite cases.
Examples for exponential overshooting can be found in the following test_suite cases:

  * 1.4M_ms_op_mono
  * 25M_pre_ms_to_core_collapse
  * 5M_cepheid_blue_loop/inlist_cepheid_blue_loop
  * 7M_prems_to_AGB/inlist_7M_prems_to_AGB
  * accretion_with_diffusion
  * agb
  * axion_cooling
  * black_hole
  * c13_pocket
  * cburn_inward
  * envelope_inflation
  * example_ccsn_IIp
  * example_make_pre_ccsn
  * gyre_in_mesa_rsg
  * high_mass
  * high_z
  * hot_cool_wind
  * magnetic_braking
  * make_co_wd
  * make_metals
  * ppisn
  * pre_zahb
  * radiative_levitation

Examples for step overshooting can be found in the following test_suite cases:

  * high_rot_darkening
  * relax_composition_j_entropy


Version number
~~~~~~~~~~~~~~

The MESA ``version_number`` is now represented as a string internally
and in the headers of history/profile output.  User scripts that
assume this is an integer may need to be revised.

``other_wind`` hook
~~~~~~~~~~~~~~~~~~~

The interface of the ``other_wind`` hook changed from ::

    subroutine other_wind_interface(id, Lsurf, Msurf, Rsurf, Tsurf, w, ierr)
       use const_def, only: dp
       integer, intent(in) :: id
       real(dp), intent(in) :: Lsurf, Msurf, Rsurf, Tsurf ! surface values (cgs)
       real(dp), intent(out) :: w ! wind in units of Msun/year (value is >= 0)
       integer, intent(out) :: ierr
    end subroutine other_wind_interface

to ::

    subroutine other_wind_interface(id, Lsurf, Msurf, Rsurf, Tsurf, X, Y, Z, w, ierr)
       use const_def, only: dp
       integer, intent(in) :: id
       real(dp), intent(in) :: Lsurf, Msurf, Rsurf, Tsurf, X, Y, Z ! surface values (cgs)
       real(dp), intent(out) :: w ! wind in units of Msun/year (value is >= 0)
       integer, intent(out) :: ierr
    end subroutine other_wind_interface

Existing user routines will need to be updated.


Removal of ``id_extra`` from ``run_star_extras.f``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Most routines in ``run_star_extras.f`` had an argument ``id_extra``.
This argument generally did not do anything and so has been removed.
Existing user routines will need to be updated.

A simple way to migrate from routines written for previous versions of
MESA is to find and replace the string ", id_extra" with the empty
string in ``run_star_extras.f``.

Change of ``extras_startup`` from function to subroutine
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The interface of ``extras_startup`` changed from ``integer function`` to subroutine.  The current empty version of this routine is::

      subroutine extras_startup(id, restart, ierr)
         integer, intent(in) :: id
         logical, intent(in) :: restart
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
      end subroutine extras_startup

Existing user routines will need to be updated to reflect this new interface.


Hooks for extra header items
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The interface of the routines

+ ``how_many_extra_history_header_items``
+ ``data_for_extra_history_header_items``
+ ``how_many_extra_profile_header_items``
+ ``data_for_extra_profile_header_items``

has changed.  If these routines are included in your
``run_star_extras.f`` (even if they have not been customized), you
will need to update them.  You should replace the old versions with::

      integer function how_many_extra_history_header_items(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_history_header_items = 0
      end function how_many_extra_history_header_items


      subroutine data_for_extra_history_header_items(id, n, names, vals, ierr)
         integer, intent(in) :: id, n
         character (len=maxlen_history_column_name) :: names(n)
         real(dp) :: vals(n)
         type(star_info), pointer :: s
         integer, intent(out) :: ierr
         ierr = 0
         call star_ptr(id,s,ierr)
         if(ierr/=0) return

         ! here is an example for adding an extra history header item
         ! also set how_many_extra_history_header_items
         ! names(1) = 'mixing_length_alpha'
         ! vals(1) = s% mixing_length_alpha

      end subroutine data_for_extra_history_header_items


      integer function how_many_extra_profile_header_items(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_profile_header_items = 0
      end function how_many_extra_profile_header_items


      subroutine data_for_extra_profile_header_items(id, n, names, vals, ierr)
         integer, intent(in) :: id, n
         character (len=maxlen_profile_column_name) :: names(n)
         real(dp) :: vals(n)
         type(star_info), pointer :: s
         integer, intent(out) :: ierr
         ierr = 0
         call star_ptr(id,s,ierr)
         if(ierr/=0) return

         ! here is an example for adding an extra profile header item
         ! also set how_many_extra_profile_header_items
         ! names(1) = 'mixing_length_alpha'
         ! vals(1) = s% mixing_length_alpha

      end subroutine data_for_extra_profile_header_items


Removal of inlist_massive_defaults
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The file inlist_massive_defaults has been removed from ``star``.
Copies of the inlist can now be found in the following test cases:

  * 25M_pre_ms_to_core_collapse
  * 25M_z2m2_high_rotation
  * adjust_net
  * black_hole
  * envelope_inflation
  * example_ccsn_IIp
  * example_make_pre_ccsn
  * magnetic_braking
  * split_burn_20M_si_burn_qp
  * split_burn_big_net_30M
  * split_burn_big_net_30M_logT_9.8

.. _Other changes r12278:

Other changes
-------------
  
* The routines ``{alloc,move,store,unpack}_extra_info`` were removed
  from ``standard_run_star_extras.inc``.  (These routines were used to
  store/retrieve information from photos.)  If you have existing
  ``run_star_extras`` code that includes these routines, it will
  continue to function.  However, in new ``run_star_extras`` code, the
  recommended way to store/retrieve data is using the
  ``other_photo_read`` and ``other_photo_write`` hooks.  Examples can
  be found in the :ref:`conductive_flame` and `brown_dwarf` test
  suite cases.

* The controls ``xtra_coef_os_*`` and ``xtra_dist_os_*`` which could
  be used to modify ``mesh_delta_coeff`` in overshooting regions have
  been removed.  The same functionality is available using the
  ``other_mesh_delta_coeff_factor`` and an example implementation is
  given in the ``agb`` test suite case.

* The output-related control ``alpha_bdy_core_overshooting`` and
  related history options ``core_overshoot_{Hp,f,f0,hstep,r0}`` and
  ``{mass,radius}_bdy_core_overshooting`` have been removed.

* The ``star_data`` module was split out of the ``star`` module.  The
  source file describing the contents of the ``star_info`` data
  structure is now located at ``star_data/public/star_data.inc``.

* If not provided with an argument, the ``./re`` script will now
  restart from the most recent photo (determined by filesystem
  modification time).

* Added star_control pre_ms_relax_to_start_radiative_core to existing
  star_control pre_ms_relax_num_steps to provide option for creating a
  pre-main sequence model just after the end of the fully convective period.   
  The relaxation steps from raw pre-ms model to end of fully convective are
  done using simple control setting selected for robustness.  After the
  relaxation is complete, the actual inlist parameter settings are used.
  
* Added a new hook other_accreting_state to allow the user to specify the
  specific total energy, pressure, and density of the accreting material.
  These properties are used by eps_mdot to compute the contribution of
  accretion to the energy equation. By default (when this hook is not used),
  these properties are all taken from the surface cell.


Changes in r12115
=================

This section describes changes that occurred since r11554.  The
changes were originally described by `this post
<https://lists.mesastar.org/pipermail/mesa-users/2019-September/010470.html>`__
to the MESA Users' mailing list.

Backwards incompatible changes
------------------------------

Changes to atmospheres
~~~~~~~~~~~~~~~~~~~~~~

There has been a major overhaul of the atmosphere controls and related
code.  This improves consistency between the atmosphere and interior
calculations and offers more flexibility to users.  To learn more,
please consult the user guide available :download:`here <assets/atm-user-guide.txt>`.

Changes to ``s% xtra`` variables
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The MESA star pointer provides a set of extra variables that can be used
in run_star_extras.f and are automatically saved and restored during
retries and backups. The old variables were

* ``s% xtra1``, ``s% xtra2``, ..., etc. for floats,
* ``s% ixtra1``, ``s% ixtra2``, ..., etc. for integers, and
* ``s% lxtra1``, ``s% lxtra2``, ..., etc. for logicals (booleans).

These have now been collapsed into arrays (e.g., ``s% xtra(:)``). If you use
these variables in your ``run_star_extras.f``, you will need to enclose the
variable number in brackets. E.g., ``s% xtra1`` becomes ``s% xtra(1)``,
``s% ixtra17`` becomes ``s% ixtra(17)``, etc.

The new scheme allows you to define integers with meaningful names that
can make it more obvious how an ``xtra`` variable is used. For example, if
you end up storing some integrated quantity in ``s% xtra(11)``, you could
define ``i_my_integral = 11`` and then refer to the value as
``s% xtra(i_my_integral)``.

The ``ppisn`` test suite case provides an example of this usage.

Other changes
-------------

Changes to WD ``atm`` tables
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

There are now 2 options for white dwarf atmosphere tables:

* ``WD_tau_25``: the original WD atmosphere table option for DA (H atmosphere)
  WDs; also found and fixed a bug in the header of this file that was
  causing it to use only a small portion of the actual table

* ``DB_WD_tau_25``: new table for DB (He atmosphere) WDs


Changes to header format
~~~~~~~~~~~~~~~~~~~~~~~~

The header format is now taken from the ``star_history_*_format`` and
``profile_*_format`` variables defined in ``controls.defaults``.  This
addresses the bug caused by the compiler version string exceeding the
allowed length of a header column found by some users with the MESA
SDK and running on macOS.  The default is now 40 characters but this
can be set to a larger (or smaller) value in ``&controls``.


In analogy to the routines in ``run_star_extras.f``,
``run_binary_extras.f`` now has the routines

::

    how_many_extra_binary_history_header_items
    data_for_extra_binary_history_header_items

that allow the user to add custom header items to the binary history
output.


New overshooting controls
~~~~~~~~~~~~~~~~~~~~~~~~~

We have introduced new, easier to use controls for overshooting, based
on convection-zone matching criteria.

Use ``overshoot_new = .true.`` to use the new controls.

Note that in a future release, these new controls will become the
default.  Therefore, when you start new inlists, we recommend that you
use these new controls.

In the new set of controls, for each convective boundary it is possible
to define an ``overshoot_zone_type``, ``overshoot_zone_loc`` and an
``overshoot_bdy_loc`` as well as values for the overshooting parameters.

The permitted values are the following:

* ``overshoot_scheme``: ``'exponential'``, ``'double_exponential'`` or ``'step'``
* ``overshoot_zone_type``: ``'burn_H'``, ``'burn_He'``, ``'burn_Z'``, ``'nonburn'`` or ``'any'``
* ``overshoot_zone_loc``: ``'core'``, ``'shell'`` or ``'any'``
* ``overshoot_bdy_loc``: ``'bottom'``, ``'top'`` or ``'any'``

The following controls assign values for the diffusive or step
overshooting parameters:

* ``overshoot_f``
* ``overshoot_f0``
* ``overshoot_f2``

The following example will apply exponential overshoot, with f = 0.128
and f0 = 0.100, at the bottom of non-burning convective shells; and
exponential overshoot, with f = 0.014 and f0 = 0.004, at all other
convective boundaries.

::

    overshoot_scheme(1) = 'exponential'
    overshoot_zone_type(1) = 'nonburn'
    overshoot_zone_loc(1) = 'shell'
    overshoot_bdy_loc(1) = 'bottom'
    overshoot_f(1) = 0.128
    overshoot_f0(1) = 0.100

    overshoot_scheme(2) = 'exponential'
    overshoot_zone_type(2) = 'any'
    overshoot_zone_loc(2) = 'any'
    overshoot_bdy_loc(2) = 'any'
    overshoot_f(2) = 0.014
    overshoot_f0(2) = 0.004

Other examples are illustrated in the ``gyre_in_mesa_rsg`` and
``high_mass`` test_suite cases.


Changes in r11554
=================

This section describes changes that occurred since r11532.  The
changes were originally described by `this post
<https://lists.mesastar.org/pipermail/mesa-users/2019-March/009905.html>`__
to the MESA Users' mailing list.

The release was principally made to quickly fix some memory leaks in r11532.
Several users saw long-running jobs killed due to
exhaustion of system memory.  Thanks to Avishai Gilkis for the report.

This release also sets the ``star_job`` control
``num_steps_for_garbage_collection = 1000``.  Periodically MESA will free some
memory from data structures that are no longer needed but have not been
deallocated yet.  At present, this only targets the EOS tables.
(Implemented by Rob Farmer)

The header of MESA history/profile files now includes information about the
compiler used and start date of the MESA run. (Implemented by Aaron Dotter)

::

               1           2        3                        4           5
  version_number    compiler    build         MESA_SDK_version        date
           11554  "gfortran"  "8.3.0"  "x86_64-linux-20190313"  "20190314"

Changes in r11532
=================

This section describes changes that occurred since r10398.  The
changes were originally described by `this post
<https://lists.mesastar.org/pipermail/mesa-users/2019-March/009842.html>`__
to the MESA Users' mailing list.

RSP is a new functionality in MESAstar that models the non-linear radial
stellar pulsations that characterize RR Lyrae, Cepheids, and other classes
of variable stars. See the ``rsp_*`` examples in the test suite.

We significantly enhance numerical energy conservation capabilities,
including during mass changes. For example, this enables calculations
through the He flash that conserve energy to better than 0.001%. Most test
cases now have this enabled, for instance ``1.3M_ms_high_Z``,
``25M_pre_ms_to_core_collapse``, and ``wd`` as examples.

To improve the modeling of rotating stars in MESA, we introduce a new
approach to modifying the pressure and temperature equations of stellar
structure, and a formulation of the projection effects of gravity
darkening. The latter are controlled by the ``grav_dark`` options in
``history_columns.list``;  see ``high_rot_darkening`` for an example of its use.

A new scheme for tracking convective boundaries, called Convective
Pre-Mixing (CPM), yields reliable values of the convective-core mass, and
allows the natural emergence of adiabatic semiconvection regions during
both core hydrogen- and helium-burning phases. Examples for this can be
found in the inlists provided with the mesa 5 paper.

We have updated the equation of state and nuclear reaction physics modules.

There are an increased number of warnings for when MESA goes beyond the
validity of the input physics (for instance the nuclear reactions rates
from REACLIB are ill-defined when logT>10.0). These warnings are controlled
by the ``warn_*`` options.

The definition of ``eps_nuc`` has slightly changed (see MESA V, Section 3.2) in
order to be suitable for use with the new energy equation.  If you are
running models using the ``dLdm`` form that includes ``eps_grav``, you should
consult the controls option ``include_composition_in_eps_grav`` and its
associated documentation.

A new set of tests (``gyre_in_mesa_*``) demonstrate how to call GYRE on the fly
during a MESA run.

The ``astero`` module now allows users to define model parameters
(``my_param[123]``) that will be optimised in a similar way to the standard
options (``mass``, ``Y``, ``FeH``, ``alpha``, ``f_ov``). These are defined in the subroutine
``set_my_params`` in ``run_star_extras.f`` in a similar way to how users can define
their own observables (``my_var[123]``).

The ``astero`` module now has controls ``normalize_chi2_*`` that allow the user to
decide whether or not to normalize each component of |chi^2| by the number of
terms that contributed to that component.

The format of the OP_MONO opacity table cache has changed.  If you have
used these files in a previous version of MESA then you should do:

::

   rm $MESA_OP_MONO_DATA_CACHE_FILENAME

before installing MESA.  If you use multiple MESA versions, this means that
you cannot share the cache file between old and new versions.  Therefore,
you should make sure to use a different cache file in each case. This may
be more easily accomplished using the controls option
``op_mono_data_cache_filename`` rather than the environment variable.

The version of GYRE bundled with MESA has been updated to version 5.2.

Binaries can now model “twins”, where we can skip the calculation of the
companion as its assumed to be identical to the primary. This is controlled
by the ``binary_job`` parameter ``*_model_twins_flag``.

There is a new way to treat convection in a model, via the
``convective_velocity_flag``. This adds an equation to solve the velocity of
convective motion, instead of using the value derived from MLT. This is
useful for models evolving on fast timescales and is a replacement for
``min_T_for_acceleration_limited_conv_velocity``.

Two new test cases (``hydro_Ttau_solar`` and ``hydro_Ttau_evolve``) demonstrate the
use of mixing length parameters and T(τ) relations calibrated to 3D
radiation-coupled hydrodynamics (RHD) simulations computed by
`Trampedach et al. (2014) <https://ui.adsabs.harvard.edu/abs/2014MNRAS.442..805T>`_. More details are provided in
`Mosumgaard et al. (2018) <https://ui.adsabs.harvard.edu/abs/2018MNRAS.478.5650M>`_. MESA
also includes low-temperature opacity tables that match those used in the
3D RHD simulations, which can be used by setting ``kappa_lowT_prefix =
'lowT_rt14_ag89'``.

There have been many bug fixes and performance enhancements to MESA.
Reports of bugs or suggested improvements are welcome on the mesa-users
mailing list.

A reminder to please share your inlists and run_star_extras on mesastar.org
upon publication of your science papers!


Changes in r10398
=================

This section describes changes that occurred since r11554.  The
changes were originally described by `this post
<https://lists.mesastar.org/pipermail/mesa-users/2018-March/008812.html>`_
to the MESA Users' mailing list.

Equation of State: PTEH, DT2, and ELM (Bill)
--------------------------------------------

Several new options for the mesa/eos have been added, all aiming for
more accurate partials for the Newton solver.  All of these new eos
options use bicubic spline interpolation in tables of lnPgas, lnS, and
lnE as a way to get numerically accurate 1st and 2nd partial
derivatives with respect to lnRho and lnT.  The partials are directly
calculated from the interpolating bicubic polynomials to give
numerical accuracy, but this comes at a cost in thermodynamic
consistency since the actual thermodynamic relations can only be
approximated by bicubic splines.

The new eos options are called "PTEH", "DT2", and "ELM".  The PTEH
tables are created using the approach of
`Pols, Tout, Eggleton, and Han (1995) <https://ui.adsabs.harvard.edu/abs/1995MNRAS.274..964P>`_
as implemented by `Paxton (2004) <https://ui.adsabs.harvard.edu/abs/2004PASP..116..699P/abstract>`_
in a program derived from
Eggleton’s stellar evolution code (1973).  PTEH extends the mesa/eos
coverage to lower densities than allowed by OPAL (down to 10^-18 g
cm^-3) and higher metallicity than covered (OPAL stops at Z = 0.04
while PTEH covers all Z).  When PTEH is enabled, it is used for low
densities and for high Z in cases that for lower Z would be handled
using data from OPAL/SCVH tables.  In the old MESA EOS we fell back to
HELM to provide approximate results for the cases now covered by PTEH.

The mesa/star default controls enable PTEH for both low densities and high Z.

::

      use_eosPTEH_for_low_density = .true.
      use_eosPTEH_for_high_Z = .true.
      Z_for_all_PTEH = 0.040d0
      Z_for_any_PTEH = 0.039d0

The two remaining new eos options, DT2 and ELM, provide high
resolution tables in logRho and logT for values from mesa/eos for
OPAL/SCVH values and for HELM respectively.  These cover a subset of
the standard eos domain with standard eos results for logPgas, logS,
and logE in a form suitable for bicubic spline interpolation in order
to give 1st and 2nd partials with high numerical accuracy. However,
since use of DT2 and ELM will give decreased thermodynamic consistency
that might not be compensated for by better residuals, these are both
disabled by default in mesa/star.

::

      use_eosDT2 = .false.
      use_eosELM = .false.

Opacities (Josiah, Aaron)
-------------------------

The opacity module (``kap``) underwent some internal restructuring.
The ``kap`` module now exposes only a single ``kap_get`` interface
instead of separate ``kap_get_Type1`` and ``kap_get_Type2``
subroutines.  This has two user-visible consequences.

* The control ``kappa_type2_logT_lower_bdy`` was removed.  That
  control was no longer needed, as the existing control
  ``kappa_blend_logT_lower_bdy`` now also applies to Type2 opacities.  All
  other related opacity controls (e.g., ``use_Type2_opacities``) remain
  unchanged in name and behavior.
* Previously, there were separate "other" hooks for Type1 and Type2
  opacities.  Now, there is only one hook, ``other_kap_get``.  It has the
  call signature of the previous Type2 hook, which is a super-set of
  the arguments to the Type1 hook (see ``star/other/other_kap.f90``).

In previous versions opacities where clipped to the edge values of the
tables when logR=logRho-3logT+18<-8. This has been replaced for a
blend to Compton opacities between logR=-7.5 and logR=-8.


Element Diffusion (Evan)
------------------------

Fixed a bug in the ionization treatment for diffusion in the pressure
ionization routine. This was due to a typo in the original paper that
presented the ionization scheme. Restored the missing factor of
rho^1/3 thanks to a later presentation of this same scheme (Dupuis et
al. 1992) and a note `here
<http://www1.astrophysik.uni-kiel.de/~koester/astrophysics/astrophysics.html>`_.

Added a user control (``D_mix_ignore_diffusion``) for when to ignore
element diffusion in surface or core mixing regions. Previously,
diffusion would be turned off for surface mixing regions of ANY
strength, even very weak mixing where diffusion might still be
relevant. Now this control is set to a D_mix of 10⁵ (cm²/s), so that
mixing that will obviously overwhelm diffusion (like convection) will
turn it off, but weaker mixing won't.


Gravity Darkening (Aaron)
-------------------------

Added options to include gravity darkening, in the form of projected (surface-averaged) luminosities and effective temperatures of the star viewed along the equator and pole, to the history file.  Assumes the star is an oblate spheroid; see `here <https://github.com/aarondotter/GDit>`_ for more info.

::

      grav_dark_L_polar !Lsun 
      grav_dark_Teff_polar !K
      grav_dark_L_equatorial !Lsun 
      grav_dark_Teff_equatorial !K


Isomers (Frank, Josiah, Bill)
-----------------------------

The isomers of ²⁶Al can now be added to a reaction network. To use
them, include the isomers in your network specification file. Two
examples include

::

    add_isos_and_reactions(
       neut
        h  1  1     ! hydrogen
       he  4  4    ! helium
       mg 25 25 ! magnesium
       al26-1      ! ground state
       al26-2      ! meta-stable excited state
       )

and ::

    include 'mesa_45.net'
    add_isos_and_reactions(
      al26-1
      al26-2
    )
    remove_iso(al26)

One may use either ``al26`` or ``al26-1`` and ``al26-2``. Reaction
rates for the ²⁶Al isomers with other isotopes are picked up from the
JINA reaclib file. Reaction rates for ``al26-1 <-> al26-2`` are
from `Gupta & Meyer (2001)
<http://adsabs.harvard.edu/abs/2001PhRvC..64b5805G>`_ and located in
``data/rates_data/rate_tables`` along with the new default
``rate_list.txt`` file.

User-Beware: if you want a local ``rate_tables`` directory,
http://mesa.sourceforge.net/star_job_defaults.html#rate_tables_dir
<http://mesa.sourceforge.net/star_job_defaults.html#rate_tables_dir>,
and you want the ²⁶Al isomers, then the two ``al26-1 <-> al26-2`` rate
files must be copied from their default location to your local
rate_tables directory and your local ``rate_list.txt`` modified to include
these two rates.


Installation Debugging (Rob, Josiah)
------------------------------------

There is a new command, ``$MESA_DIR/help`` which outputs system
information we need when debugging installation issues and/or MESA
crashes. ``./install`` will now also log its output to a file
``$MESA_DIR/build.log``, if you have an installation issue please include
this file when reporting an issue to mesa-users.


Miscellaneous improvements (Rob, Josiah)
----------------------------------------

You can now use the ``MESA_INLIST`` environment variable to set the
name of the main inlist file when using MESA binary.

The output cadence of MESA binary has been tweaked to that its
behavior is the same as MESA star.  (If you use the same options, you
should get output at the same steps.)

There is now a flag ``b% need_to_update_binary_history_now``, which if
set forces binary history output to occur at the current step.


Run_star_extras (Aaron)
-----------------------

Put calls to ``extra_header_items`` back into
``standard_run_star_extras`` and provided working examples of how to
call all of them.  These are useful for adding extra information to
the history and profile headers beyond what is provided by default,
such as including ``mixing_length_alpha`` in the history file header.


Building with Other Compilers
-----------------------------

MESA currently does not compile with ifort.
Other non-SDK compilers that are known to work (at the bit-for-bit level):
Gfortran 7.3.1 (fedora 27)
