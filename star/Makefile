include ../make/defaults-module.mk

# Build

MODULE_NAME := star
SRCS := \
   other/other_accreting_state.f90 \
   other/other_adjust_mdot.f90 \
   other/other_adjust_mlt_gradt_fraction.f90 \
   other/other_after_set_mixing_info.f90 \
   other/other_after_solver_setmatrix.f90 \
   other/other_alpha_mlt.f90 \
   other/other_am_mixing.f90 \
   other/other_astero_freq_corr.f90 \
   other/other_brunt.f90 \
   other/other_brunt_smoothing.f90 \
   other/other_build_initial_model.f90 \
   other/other_cgrav.f90 \
   other/other_close_gaps.f90 \
   other/other_d_mix.f90 \
   other/other_diffusion.f90 \
   other/other_diffusion_coefficients.f90 \
   other/other_diffusion_factor.f90 \
   other/other_energy.f90 \
   other/other_energy_implicit.f90 \
   other/other_eps_grav.f90 \
   other/other_eval_fp_ft.f90 \
   other/other_eval_i_rot.f90 \
   job/other_extras.f90 \
   other/other_gradr_factor.f90 \
   other/other_j_for_adjust_j_lost.f90 \
   other/other_kap.f90 \
   other/other_mesh_delta_coeff_factor.f90 \
   other/other_mesh_functions.f90 \
   other/other_mlt_results.f90 \
   other/other_momentum.f90 \
   other/other_momentum_implicit.f90 \
   other/other_net_derivs.f90 \
   other/other_net_get.f90 \
   other/other_neu.f90 \
   other/other_opacity_factor.f90 \
   other/other_overshooting_scheme.f90 \
   other/other_pgstar_plots.f90 \
   other/other_photo_read.f90 \
   other/other_photo_write.f90 \
   other/other_pressure.f90 \
   other/other_rate_get.f90 \
   other/other_remove_surface.f90 \
   other/other_rsp_build_model.f90 \
   other/other_rsp_linear_analysis.f90 \
   other/other_screening.f90 \
   other/other_set_pgstar_controls.f90 \
   other/other_solver_monitor.f90 \
   other/other_split_burn.f90 \
   other/other_surface_pt.f90 \
   other/other_timestep_limit.f90 \
   other/other_torque.f90 \
   other/other_torque_implicit.f90 \
   other/other_wind.f90 \
   private/adjust_mass.f90 \
   private/adjust_mesh.f90 \
   private/adjust_mesh_split_merge.f90 \
   private/adjust_mesh_support.f90 \
   private/adjust_net.f90 \
   private/adjust_xyz.f90 \
   private/alloc.f90 \
   private/atm_support.f90 \
   private/auto_diff_support.f90 \
   private/brunt.f90 \
   private/conv_premix.f90 \
   private/create_initial_model.f90 \
   private/ctrls_io.f90 \
   private/diffusion.f90 \
   private/diffusion_procs.f90 \
   private/diffusion_support.f90 \
   private/do_one_utils.f90 \
   private/element_diffusion.f90 \
   private/eos_support.f90 \
   private/eps_grav.f90 \
   private/eps_mdot.f90 \
   private/evolve.f90 \
   private/evolve_support.f90 \
   private/gravity_darkening.f90 \
   private/history.f90 \
   private/history_specs.f90 \
   private/hydro_alpha_rti_eqns.f90 \
   private/hydro_chem_eqns.f90 \
   private/hydro_energy.f90 \
   private/hydro_eqns.f90 \
   private/hydro_momentum.f90 \
   private/hydro_riemann.f90 \
   private/hydro_rotation.f90 \
   private/hydro_rsp2.f90 \
   private/hydro_rsp2_support.f90 \
   private/hydro_temperature.f90 \
   private/hydro_vars.f90 \
   private/init.f90 \
   private/init_model.f90 \
   private/ionization_potentials.f90 \
   private/kap_support.f90 \
   private/magnetic_diffusion.f90 \
   private/mass_utils.f90 \
   private/mesh_adjust.f90 \
   private/mesh_functions.f90 \
   private/mesh_plan.f90 \
   private/micro.f90 \
   private/mix_info.f90 \
   private/mod_typical_charge.f90 \
   private/net.f90 \
   private/neu.f90 \
   private/overshoot.f90 \
   private/overshoot_exp.f90 \
   private/overshoot_step.f90 \
   private/overshoot_utils.f90 \
   private/paquette_coeffs.f90 \
   private/phase_separation.f90 \
   private/photo_in.f90 \
   private/photo_out.f90 \
   private/pre_ms_model.f90 \
   private/predictive_mix.f90 \
   private/profile.f90 \
   private/profile_getval.f90 \
   private/pulse.f90 \
   private/pulse_cafein.f90 \
   private/pulse_fgong.f90 \
   private/pulse_gr1d.f90 \
   private/pulse_gsm.f90 \
   private/pulse_gyre.f90 \
   private/pulse_osc.f90 \
   private/pulse_saio.f90 \
   private/pulse_utils.f90 \
   private/read_model.f90 \
   private/relax.f90 \
   private/remove_shells.f90 \
   private/report.f90 \
   private/rotation_mix_info.f90 \
   private/rsp.f90 \
   private/rsp_build.f90 \
   private/rsp_def.f90 \
   private/rsp_eval_eos_and_kap.f90 \
   private/rsp_lina.f90 \
   private/rsp_relax_env.f90 \
   private/rsp_step.f90 \
   job/run_star_support.f90 \
   private/set_flags.f90 \
   private/solve_omega_mix.f90 \
   private/solver_support.f90 \
   private/star_bcyclic.f90 \
   private/star_history_def.f90 \
   private/star_job_ctrls_io.f90 \
   private/star_private_def.f90 \
   private/star_profile_def.f90 \
   private/star_solver.f90 \
   private/star_utils.f90 \
   private/starspots.f90 \
   private/struct_burn_mix.f90 \
   private/timestep.f90 \
   private/turb_info.f90 \
   private/turb_support.f90 \
   private/winds.f90 \
   private/write_model.f90 \
   public/star_def.f90 \
   public/star_lib.f90 \
   job/run_star.f90

ifeq ($(WITH_PGSTAR),yes)
SRCS += \
   other/pgstar_decorator.f90 \
   private/pgstar_ctrls_io.f90 \
   private/pgstar_support.f90 \
   private/pgstar_hist_track.f90 \
   private/pgstar_kipp.f90 \
   private/pgstar_dpg_dnu.f90 \
   private/pgstar_trho.f90 \
   private/pgstar_tmaxrho.f90 \
   private/pgstar_hr.f90 \
   private/pgstar_r_l.f90 \
   private/pgstar_r_teff.f90 \
   private/pgstar_l_teff.f90 \
   private/pgstar_l_r.f90 \
   private/pgstar_l_v.f90 \
   private/pgstar_logl_teff.f90 \
   private/pgstar_logl_r.f90 \
   private/pgstar_logl_v.f90 \
   private/pgstar_logg_teff.f90 \
   private/pgstar_logg_logt.f90 \
   private/pgstar_history_panels.f90 \
   private/pgstar_trho_profile.f90 \
   private/pgstar_mode_prop.f90 \
   private/pgstar_power.f90 \
   private/pgstar_production.f90 \
   private/pgstar_abundance.f90 \
   private/pgstar_dynamo.f90 \
   private/pgstar_mixing_ds.f90 \
   private/pgstar_summary_burn.f90 \
   private/pgstar_summary_profile.f90 \
   private/pgstar_summary_history.f90 \
   private/pgstar_summary.f90 \
   private/pgstar_color_magnitude.f90 \
   private/pgstar_profile_panels.f90 \
   private/pgstar_grid.f90 \
   private/pgstar_network.f90 \
   other/sample_pgstar_plot.f90 \
   private/pgstar_full.f90
EXTERNAL_DEPENDS_ON := pgplot
else
SRCS += private/pgstar_stub.f90
EXTERNAL_DEPENDS_ON :=
endif

SRCS_CHECK = test/src/run.f90 test/src/run_star_extras.f90
INTERNAL_DEPENDS_ON := const utils star_data ionization atm turb auto_diff
BINTYPE := static-lib
INCLUDE_DIRS := -Iprivate -Ipublic -Idefaults
CHECK_DIFF_PROG := sh -c "diff -b <(grep 'stop because star_age >= max_age' \$$1) \$$2" diff

# Testing

CHECK_RESULTS_GOLDEN := test/test_output

# Install

MODULES := star_lib.mod star_def.mod run_star_support.mod other_extras.mod pulse.mod \
           pgstar_abundance.mod \
           pgstar_color_magnitude.mod \
           pgstar_ctrls_io.mod \
           pgstar_decorator.mod \
           pgstar_dpg_dnu.mod \
           pgstar_dynamo.mod \
           pgstar_grid.mod \
           pgstar_history_panels.mod \
           pgstar_hist_track.mod \
           pgstar_hr.mod \
           pgstar_kipp.mod \
           pgstar_logg_logt.mod \
           pgstar_logg_teff.mod \
           pgstar_logl_r.mod \
           pgstar_logl_teff.mod \
           pgstar_logl_v.mod \
           pgstar_l_r.mod \
           pgstar_l_teff.mod \
           pgstar_l_v.mod \
           pgstar_mixing_ds.mod \
           pgstar.mod \
           pgstar_mode_prop.mod \
           pgstar_network.mod \
           pgstar_power.mod \
           pgstar_production.mod \
           pgstar_profile_panels.mod \
           pgstar_r_l.mod \
           pgstar_r_teff.mod \
           pgstar_summary_burn.mod \
           pgstar_summary_history.mod \
           pgstar_summary.mod \
           pgstar_summary_profile.mod \
           pgstar_support.mod \
           pgstar_tmaxrho.mod \
           pgstar_trho.mod \
           pgstar_trho_profile.mod

include $(MAKE_DIR)/Makefile
