.. _list-otherhooks:

eos hooks
=========

For details of the eos hooks, see ``eos/other``.    Code location
indicates the file in ``eos/private`` where the hook is called.

===============================  =============
Hook name                        Code Location
===============================  =============
other_eos_frac                   eosdt_eval
other_eos_component              eosdt_eval
other_eos_results                eosdt_eval
===============================  =============


star hooks
==========

There are many hooks in star (see ``star/other``).  Code location
indicates the file in ``star/private`` where the hook is called.

Equations
---------
===============================  =============
Hook name                        Code Location
===============================  =============
other_energy                     evolve
other_energy_implicit            hydro_eqns
other_momentum                   hydro_vars
other_momentum_implicit          hydro_vars
other_eps_grav                   eps_grav
other_gradr_factor               hydro_vars
other_opacity_factor             micro
other_cgrav                      hydro_vars
===============================  =============
   
atm
---
===============================  =============
Hook name                        Code Location
===============================  =============
other_surface_pt                 hydro_vars
===============================  =============

kap
---
===============================  =============
Hook name                        Code Location
===============================  =============
other_kap_get                    kap_support
other_kap_get_op_mono            kap_support
===============================  =============

net
---
===============================  =============
Hook name                        Code Location
===============================  =============
other_net_get                    net
set_rate_factors                 net
set_which_rates                  net
===============================  =============

neu
---
===============================  =============
Hook name                        Code Location
===============================  =============
other_neu                        neu
===============================  =============

diffusion
---------
===============================  =============
Hook name                        Code Location
===============================  =============
other_diffusion                  element_diffusion
other_diffusion_factor           element_diffusion
other_diffusion_coefficients     diffusion_support
===============================  =============

rotation
--------
===============================  =============
Hook name                        Code Location
===============================  =============
other_torque                     solve_omega_mix
other_torque_implicit            solve_omega_mix
binary_other_torque_implicit     solve_omega_mix
other_eval_fp_ft                 hydro_rotation
other_am_mixing                  mix_info
other_j_for_adjust_j_lost        adjust_mass
===============================  =============

composition mixing
------------------
===============================  =============
Hook name                        Code Location
===============================  =============
other_alpha_mlt                  hydro_vars
other_mlt                        mlt_info
other_d_mix                      mix_info
other_adjust_mlt_gradt_fraction  hydro_vars
other_overshooting_scheme        overshoot
other_split_mix                  solve_mix
other_after_set_mixing_info      mix_info
===============================  =============

winds and mass change
---------------------
===============================  =============
Hook name                        Code Location
===============================  =============
other_wind                       winds
other_adjust_mdot                evolve
===============================  =============

brunt and astero
----------------
===============================  =============
Hook name                        Code Location
===============================  =============
other_brunt                      brunt
other_brunt_smoothing            brunt
other_astero_freq_corr           astero_support
===============================  =============

pgstar
------
===============================  =============
Hook name                        Code Location
===============================  =============
other_pgstar_plots_info          pgstar_full
pgstar_decorator                 pgstar_support
===============================  =============

mesh
----
===============================  =============
Hook name                        Code Location
===============================  =============
how_many_other_mesh_fcns         mesh_functions
other_mesh_fcn_data              mesh_functions
other_mesh_delta_coeff_factor    adjust_mesh_support
===============================  =============

timesteps
---------
===============================  =============
Hook name                        Code Location
===============================  =============
other_timestep_limit             timestep
===============================  =============

rsp
---
===============================  =============
Hook name                        Code Location
===============================  =============
other_rsp_build_model            rsp
other_rsp_linear_analysis        rsp
===============================  =============

photos
------
===============================  =============
Hook name                        Code Location
===============================  =============
other_photo_read                 photo_in
other_photo_write                photo_out
===============================  =============

logs
----
=====================================  =============
Hook name                              Code Location
=====================================  =============
how_many_extra_history_columns         history
data_for_extra_history_columns         
how_many_extra_profile_columns         profile
data_for_extra_profile_columns
how_many_extra_history_header_items
data_for_extra_history_header_items
how_many_extra_profile_header_items
data_for_extra_profile_header_items
data_for_extra_binary_history_columns
=====================================  =============

initial model
-------------
===============================  =============
Hook name                        Code Location
===============================  =============
other_build_initial_model        create_initial_model
===============================  =============

relax
-----
===============================  =============
Hook name                        Code Location
===============================  =============
finished_relax                   relax
===============================  =============

solver
------
===============================  =============
Hook name                        Code Location
===============================  =============
other_after_enter_setmatrix      hydro_mtx		
other_after_struct_burn_mix	 struct_burn_mix
other_before_struct_burn_mix	 struct_burn_mix
other_solver_monitor		 star_solver
other_new_generation             evolve_support
other_set_current_to_old         evolve_support
===============================  =============

job extras
----------
===============================  =============
Hook name                        Code Location
===============================  =============
extras_startup                   run_star_support
extras_controls                  run_star_support
extras_check_model               run_star_support
extras_finish_step               run_star_support
extras_after_evolve              run_star_support
===============================  =============
