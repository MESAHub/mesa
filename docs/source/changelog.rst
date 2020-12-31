*********
Changelog
*********

Changes in dev
==============

Backwards-incompatible changes
------------------------------




Module-level changes
--------------------

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
         integer, intent(in) :: handle ! from alloc_kap_handle
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


neu
~~~

The call signature of other_neu has changed. You no longer need to pass in z2bar



Changes in r15140
=================

Backwards-incompatible changes
------------------------------

Addition of eos and kap namelists
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The options associated with the ``eos`` and ``kap`` modules have been
moved into their own namelists.  (That is, there now exist ``&eos``
and ``&kap`` at the same level as ``&star_job`` and ``&controls``.)
User inlists will need to be updated.  See :ref:`Module-level changes`
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
controls is described :ref:`in the developer documentation <Diagnosing Solver Struggles>`.



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
applied <Overview of eos module>` to understand if MESA provides
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
The :ref:`module controls <eos module controls>` and their default
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
<Overview of eos module>`.


gyre
~~~~

GYRE has been upgraded to version 6.0.  See the `GYRE Documentation
<https://gyre.readthedocs.io/en/latest/index.html>`__ for information
about this release.

kap
~~~

Opacity-related options have been moved into their own ``kap`` namelist.
The :ref:`module controls <kap module controls>` and their default
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
combined, see the :ref:`new overview of the kap module <Overview of
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



Other changes
-------------
  
* The routines ``{alloc,move,store,unpack}_extra_info`` were removed
  from ``standard_run_star_extras.inc``.  (These routines were used to
  store/retrieve information from photos.)  If you have existing
  ``run_star_extras`` code that includes these routines, it will
  continue to function.  However, in new ``run_star_extras`` code, the
  recommended way to store/retrieve data is using the
  ``other_photo_read`` and ``other_photo_write`` hooks.  Examples can
  be found in the :ref:`conductive_flame` and :ref:`brown_dwarf` test
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
