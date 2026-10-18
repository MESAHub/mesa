.. dev_TDC_Cepheid_Pulsation:

**************************
dev_TDC_Cepheid_Pulsation
**************************

This example loads a saved stellar model, removes the core, constructs a TDC
pulsation mesh, and applies a radial velocity perturbation from Star LNA. The
nonlinear hydrodynamic calculation then measures the period, growth rate,
radius variation, luminosity variation, and surface velocity of each cycle.
GYRE is not used by this example.

The supplied ``standard_he_dep.mod`` is the 6 M core helium burning model from
``dev_TDC_Cepheid_6M``. The supplied ``model_00000800.mod`` is the 80 M
star-top comparison model. Its commented load option is next to the active 6 M
option in ``inlist_pulses``. When selecting the 80 M model, also select the
commented GS98 opacity values, ``Zbase = 0.02d0``, and the 5 MK center-removal
temperature.

The one-time TDC envelope remesh constructs 190 zones. This retains 150
mass-grid zones, reserves 20 zones next to the inner boundary, and adds 20
temperature-gradient refinement zones. The reserved inner cells decrease in
mass to ``max_center_cell_dq = 1d-7`` at the inner boundary.
``steps_before_remesh_for_TDC_pulsations`` selects its timing. A value of
``-1`` disables the remesh, ``0`` remeshes before the first step, and a
positive value remeshes after that many accepted steps. The supplied value is
``0``. A photo written before a deferred remesh retains the original target.
The remesh is not repeated by a photo written after the target.

The supplied calculation enables standard mesh adjustment. The
``restore_mesh_on_retry`` and ``num_steps_to_hold_mesh_after_retry`` controls
are listed but disabled, so changing ``okay_to_remesh`` on a restart takes
effect immediately. Enabling retry restoration can repeatedly suppress mesh
changes in a model that retries when the hold expires.

Hydrodynamic Variable
=====================

The ``&star_job`` namelist explicitly sets both velocity flags. The supplied
configuration uses ``v_flag``::

   new_v_flag = .true.
   new_u_flag = .false.

Swap these values to use ``u_flag``. Velocity time centering and artificial
viscosity are used only by the ``v_flag`` equations. Star LNA requires the
turbulent-pressure perturbation for the supplied TDC background when using
``u_flag``. The final ``u_flag`` section contains the PPISN metric split/merge
AMR scheme, including pressure reconstruction across child cells after a
split. The supplied calculation enables standard mesh adjustment with
``okay_to_remesh = .true.`` and ``use_split_merge_amr = .false.``. Set
``use_split_merge_amr = .true.`` to use the PPISN split/merge scheme instead.

Star LNA
========

Star LNA runs at model 200, writes 15 radial modes to ``LNA``, and uses the
fundamental mode for a 5 km/s surface velocity perturbation. The calculation
includes perturbations to convective flux, turbulent pressure, eddy viscosity,
and turbulent energy. ``star_LNA_stop_after_run = .false.`` continues into the
nonlinear calculation after the modes are written.

The TDC section sets ``harmonic_dissipation_length_beta = 1``. The mixing and
dissipation length is therefore

.. math::

   \frac{1}{\Lambda}
   = \frac{1}{\alpha_{\rm MLT}H_P} + \frac{1}{\beta r}.

Timestep Resolution
===================

``x_integer_ctrl(8)`` sets the requested number of timesteps per period. Once
the model number exceeds ``x_ctrl(13)``, if this value is ``N``, the additional
timestep limit before a measured period is

.. math::

   \Delta t_{\max} = \frac{\tau_{\rm dyn}}{N}.

After the first measured period, the limit becomes

.. math::

   \Delta t_{\max} = \frac{P}{N}.

The fallback ceiling from ``x_ctrl(18)`` is also applied, and the smaller value
is used. The supplied ``1d6`` second ceiling normally leaves the dynamical or
period resolution in control. Set ``x_ctrl(18) <= 0`` to disable the fallback,
or set ``x_integer_ctrl(8) <= 0`` to disable the period-based limit.

Temperature Gradient
====================

The default calculation uses the standard MESA temperature-gradient equation.
Set ``use_dPrad_dm_form_of_T_gradient_eqn = .true.`` to use the
``dPrad/dm`` form. The associated opacity floor and flux-limiting controls are
listed in the temperature-gradient section. ``convergence_ignore_equL_residuals``
is included there as a commented option, and
``make_gradr_sticky_in_solver_iters`` controls the MLT radiative-gradient
branch during Newton iterations. With ``x_logical_ctrl(25) = .true.``,
``run_star_extras`` ignores the residual for the first ten models and includes
it afterward. Set the control false before enabling the native convergence
option.

Surface Boundary
================

The momentum outer boundary condition uses
``floor_momentum_outer_BC_at_Prad = .true.``. It applies

.. math::

   P_{\rm bc} \leftarrow \max\left(P_{\rm bc},\frac{aT_{\rm bc}^4}{3}\right)

before evaluating the momentum boundary condition. This prevents the
atmosphere from imposing a negative gas pressure. It is separate from
``min_kap_for_dPrad_dm_eqn``, which floors the face opacity in the
temperature-gradient equation.

Last-Updated: 2026-08-27 by Ebraheem Farag
