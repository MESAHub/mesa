! ***********************************************************************
!
! Copyright (C) 2010 The Mesa Team
!
! MESA is free software; you can use it and/or modify
! it under the combined terms and restrictions of the MESA MANIFESTO
! and the GNU General Library Public License as published
! by the Free Software Foundation; either version 2 of the License,
! or (at your option) any later version.
!
! You should have received a copy of the MESA MANIFESTO along with
! this software; if not, it is available at the mesa website:
! http://mesa.sourceforge.net/
!
! MESA is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
! See the GNU Library General Public License for more details.
!
! You should have received a copy of the GNU Library General Public License
! along with this software; if not, write to the Free Software
! Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
!
! ***********************************************************************

 module ctrls_io

 use const_def
 use star_private_def

 implicit none

 include 'star_controls.inc'
 include 'star_controls_dev.inc'

 logical :: read_extra_controls_inlist1
 character (len=strlen) :: extra_controls_inlist1_name

 logical :: read_extra_controls_inlist2
 character (len=strlen) :: extra_controls_inlist2_name

 logical :: read_extra_controls_inlist3
 character (len=strlen) :: extra_controls_inlist3_name

 logical :: read_extra_controls_inlist4
 character (len=strlen) :: extra_controls_inlist4_name

 logical :: read_extra_controls_inlist5
 character (len=strlen) :: extra_controls_inlist5_name

 logical :: save_controls_namelist
 character (len=strlen) :: controls_namelist_name

 namelist /controls/ &
 
    ! where to start
    initial_mass, initial_z, initial_y, initial_he3, &
    
    ! definition of core boundaries
    he_core_boundary_h1_fraction, co_core_boundary_he4_fraction, one_core_boundary_he4_c12_fraction, &
    fe_core_boundary_si28_fraction, neutron_rich_core_boundary_Ye_max, min_boundary_fraction, &
    
    ! when to stop
    max_model_number, relax_max_number_retries, max_number_retries, max_age, max_age_in_seconds, max_age_in_days, &
    num_adjusted_dt_steps_before_max_age, dt_years_for_steps_before_max_age, max_abs_rel_run_E_err, &
    reduction_factor_for_max_timestep, when_to_stop_rtol, when_to_stop_atol, &
    gamma_center_limit, eta_center_limit, log_center_temp_limit, log_max_temp_upper_limit, &
    log_max_temp_lower_limit , log_center_temp_lower_limit, log_center_density_limit, &
    log_center_density_lower_limit, center_entropy_limit, center_entropy_lower_limit, &
    max_entropy_limit, max_entropy_lower_limit, min_timestep_limit, non_fe_core_rebound_limit, &
    fe_core_infall_limit, center_Ye_lower_limit, center_R_lower_limit, non_fe_core_infall_limit, &
    v_div_csound_surf_limit, v_div_csound_max_limit, Lnuc_div_L_upper_limit, Lnuc_div_L_lower_limit,&
    v_surf_div_v_kh_upper_limit, v_surf_div_v_kh_lower_limit, v_surf_div_v_esc_limit, v_surf_kms_limit, &
    stop_near_zams, Lnuc_div_L_zams_limit, Pgas_div_P_limit, Pgas_div_P_limit_max_q, gamma1_limit, gamma1_limit_max_q, &
    stop_at_phase_PreMS, stop_at_phase_ZAMS, stop_at_phase_IAMS, stop_at_phase_TAMS, gamma1_limit_max_v_div_vesc, &
    stop_at_phase_He_Burn, stop_at_phase_ZACHeB, stop_at_phase_TACHeB, &
    stop_at_phase_TP_AGB, stop_at_phase_C_Burn, stop_at_phase_Ne_Burn, &
    stop_at_phase_O_Burn, stop_at_phase_Si_Burn, stop_at_phase_WDCS, &
    peak_burn_vconv_div_cs_limit, omega_div_omega_crit_limit, delta_nu_lower_limit, &
    delta_nu_upper_limit, delta_Pg_lower_limit, delta_Pg_upper_limit, shock_mass_upper_limit, &
    mach1_mass_upper_limit, stop_when_reach_this_cumulative_extra_heating, &
    xa_central_lower_limit_species, xa_central_lower_limit, xa_central_upper_limit_species, xa_central_upper_limit, &
    xa_surface_lower_limit_species, xa_surface_lower_limit, xa_surface_upper_limit_species, xa_surface_upper_limit, &
    xa_average_lower_limit_species, xa_average_lower_limit, xa_average_upper_limit_species, xa_average_upper_limit, &
    star_species_mass_min_limit, star_species_mass_min_limit_iso, star_species_mass_max_limit, star_species_mass_max_limit_iso, &
    xmstar_min_limit, xmstar_max_limit, envelope_mass_limit, envelope_fraction_left_limit, &
    he_core_mass_limit, co_core_mass_limit, one_core_mass_limit, &
    fe_core_mass_limit, neutron_rich_core_mass_limit, HB_limit, star_mass_min_limit, star_mass_max_limit, &
    he_layer_mass_lower_limit, abs_diff_lg_LH_lg_Ls_limit, Teff_upper_limit, Teff_lower_limit, &
    photosphere_m_upper_limit, photosphere_m_lower_limit, photosphere_m_sub_M_center_limit, &
    photosphere_r_upper_limit, photosphere_r_lower_limit, log_Teff_upper_limit, log_Teff_lower_limit, &
    log_Tsurf_upper_limit, log_Tsurf_lower_limit, log_Rsurf_upper_limit, log_Rsurf_lower_limit, &
    log_Psurf_upper_limit, log_Psurf_lower_limit, remnant_mass_min_limit, ejecta_mass_max_limit, &
    log_Dsurf_upper_limit, log_Dsurf_lower_limit, log_L_upper_limit, log_L_lower_limit, &
    log_g_upper_limit, log_g_lower_limit, power_nuc_burn_upper_limit, power_h_burn_upper_limit, &
    power_he_burn_upper_limit, power_z_burn_upper_limit, power_nuc_burn_lower_limit, &
    power_h_burn_lower_limit, power_he_burn_lower_limit, power_z_burn_lower_limit, &
    
    ! max timesteps
    max_timestep, max_years_for_timestep, &
    hi_T_max_years_for_timestep, max_timestep_hi_T_limit, &
    
    ! output of "snapshots" for restarts
    photo_interval, photo_digits, photo_directory, &
    
    ! output of logs and profiles
    do_history_file, history_interval, write_header_frequency, terminal_interval, &
    terminal_show_age_units, terminal_show_timestep_units, terminal_show_log_dt, terminal_show_log_age, &
    num_trace_history_values, trace_history_value_name, write_profiles_flag, profile_interval, &
    priority_profile_interval, log_directory, star_history_name, star_history_header_name, &
    star_history_dbl_format, star_history_int_format, star_history_txt_format, extra_terminal_output_file, &
    profiles_index_name, profile_data_prefix, profile_data_suffix, profile_data_header_suffix, &
    profile_int_format, profile_txt_format, profile_dbl_format, profile_header_include_sys_details, &
    profile_model, max_num_profile_models, max_num_profile_zones, &
    write_controls_info_with_profile, controls_data_prefix, controls_data_suffix, &
    write_pulse_data_with_profile, pulse_data_format, add_atmosphere_to_pulse_data, &
    add_center_point_to_pulse_data, keep_surface_point_for_pulse_data, add_double_points_to_pulse_data, &
    interpolate_rho_for_pulse_data, threshold_grad_mu_for_double_point, max_number_of_double_points,&
    fgong_header, fgong_ivers, &
    max_num_gyre_points, format_for_OSC_data, &
    fgong_zero_A_inside_r, use_other_export_pulse_data, use_other_get_pulse_data, use_other_edit_pulse_data, &
    write_model_with_profile, model_data_prefix, model_data_suffix, &
    mixing_D_limit_for_log, trace_mass_location, min_tau_for_max_abs_v_location, &
    min_q_for_inner_mach1_location, max_q_for_outer_mach1_location, &
    mass_depth_for_L_surf, conv_core_gap_dq_limit, &
    max_dt_div_tau_conv_for_TDC, max_dt_years_for_TDC, max_X_for_gradT_eqn, &
    alpha_TDC_DAMP, alpha_TDC_DAMPR, alpha_TDC_PtdVdt, &
    
    ! burn zone eps definitions for use in logs and profiles
    burn_min1, burn_min2, &
    max_conv_vel_div_csound_maxq, width_for_limit_conv_vel, max_q_for_limit_conv_vel, &
    max_mass_in_gm_for_limit_conv_vel, max_r_in_cm_for_limit_conv_vel, &
    
    ! for reported surface/center abundances
    surface_avg_abundance_dq, center_avg_value_dq, &
    
    ! mixing parameters 
    min_convective_gap, min_thermohaline_gap, min_semiconvection_gap, min_thermohaline_dropout, &
    max_dropout_gradL_sub_grada, remove_embedded_semiconvection, recalc_mix_info_after_evolve, remove_mixing_glitches, &
    okay_to_remove_mixing_singleton, prune_bad_cz_min_Hp_height, prune_bad_cz_min_log_eps_nuc, &
    redo_conv_for_dr_lt_mixing_length, alpha_semiconvection, okay_to_reduce_gradT_excess, &
    semiconvection_option, use_Ledoux_criterion, D_mix_zero_region_bottom_q, &
    num_cells_for_smooth_gradL_composition_term, threshold_for_smooth_gradL_composition_term, clip_D_limit, &
   gradT_excess_f1, gradT_excess_f2, gradT_excess_age_fraction, gradT_excess_max_change, gradT_excess_lambda1, &
   gradT_excess_beta1, gradT_excess_lambda2, gradT_excess_beta2, gradT_excess_dlambda, gradT_excess_dbeta, gradT_excess_max_center_h1, &
   gradT_excess_min_center_he4, gradT_excess_max_logT, gradT_excess_min_log_tau_full_on, gradT_excess_max_log_tau_full_off, &
    use_superad_reduction, superad_reduction_gamma_limit, superad_reduction_gamma_limit_scale, D_mix_zero_region_top_q, &
    superad_reduction_gamma_inv_scale, superad_reduction_diff_grads_limit, superad_reduction_limit, &
    fix_eps_grav_transition_to_grid, make_gradr_sticky_in_solver_iters, min_logT_for_make_gradr_sticky_in_solver_iters, &
    max_logT_for_mlt, thermohaline_coeff, thermohaline_option, mixing_length_alpha, remove_small_D_limit, &
    alt_scale_height_flag, Henyey_MLT_y_param, Henyey_MLT_nu_param, no_MLT_below_shock, mlt_make_surface_no_mixing, &
    MLT_option, mlt_use_rotation_correction, mlt_Pturb_factor, do_normalize_dqs_as_part_of_set_qs, &
    max_Y_for_burn_z_mix_region, max_X_for_burn_he_mix_region, &
    limit_overshoot_Hp_using_size_of_convection_zone, RSP_min_tau_for_turbulent_flux, &
    predictive_mix, predictive_superad_thresh, predictive_avoid_reversal, predictive_limit_ingestion,&
    predictive_ingestion_factor, predictive_zone_type, predictive_zone_loc, predictive_bdy_loc, &
    predictive_bdy_q_min, predictive_bdy_q_max, T_mix_limit, RSP_report_undercorrections, &
    do_conv_premix, conv_premix_avoid_increase, conv_premix_time_factor, &
    conv_premix_fix_pgas, conv_premix_dump_snapshots, do_premix_heating, &
    overshoot_f, overshoot_f0, overshoot_D0, RSP_Qvisc_linear, dq_D_mix_zero_at_H_He_crossover, &
    overshoot_Delta0, overshoot_mass_full_on, overshoot_mass_full_off, dq_D_mix_zero_at_H_C_crossover, &
    overshoot_scheme, overshoot_zone_type, overshoot_zone_loc, RSP_Qvisc_quadratic, &
    overshoot_bdy_loc, overshoot_D_min, overshoot_brunt_B_max, mlt_gradT_fraction, max_conv_vel_div_csound, &
    max_v_for_convection, max_q_for_convection_with_hydro_on, alpha_RTI_src_max_q, &
    max_v_div_cs_for_convection, max_abs_du_div_cs_for_convection, RSP_max_dt, RSP_relax_dm_tolerance, &
    calculate_Brunt_B, calculate_Brunt_N2, brunt_N2_coefficient, num_cells_for_smooth_brunt_B, &
    threshold_for_smooth_brunt_B, min_magnitude_brunt_B, RSP_max_dt_times_min_rad_diff_time, &
    min_overshoot_q, overshoot_alpha, RSP_target_steps_per_cycle, &
    RSP_max_num_periods, RSP_min_max_R_for_periods, RSP_min_deltaR_for_periods, &
    RSP_min_PERIOD_div_PERIODLIN, RSP_report_limit_dt, RSP_mode_for_setting_PERIODLIN, RSP_initial_dt_factor, &
    RSP_v_div_cs_threshold_for_dt_limit, RSP_max_dt_times_min_dr_div_cs, RSP_thetae, &
    RSP_alfa, RSP_thetaq, RSP_default_PERIODLIN, &
   RSP_theta, RSP_thetat, RSP_thetau, RSP_wtr, RSP_wtc, RSP_wtt, RSP_gam, RSP_max_retries_per_step, RSP_Qvisc_linear_static, &
   RSP_alfa, RSP_alfap, RSP_alfam, RSP_alfat, RSP_alfas, RSP_alfac, RSP_alfad, RSP_gammar, RSP_nz_div_IBOTOM, &
   RSP_efl0, RSP_cq, RSP_zsh, RSP_tol_max_corr, RSP_tol_max_resid, RSP_max_iters_per_try, &
    RTI_smooth_mass, RTI_smooth_iterations, RTI_smooth_fraction, RSP_dq_1_factor, &
    alpha_RTI_diffusion_factor, dudt_RTI_diffusion_factor, dedt_RTI_diffusion_factor, alpha_RTI_src_min_v_div_cs, &
    dlnddt_RTI_diffusion_factor, composition_RTI_diffusion_factor, max_M_RTI_factors_full_on, min_M_RTI_factors_full_off, &
    alpha_RTI_src_min_v_div_cs, alpha_RTI_src_max_q, alpha_RTI_src_min_q, RSP_kick_vsurf_km_per_sec, &
    RSP_nz_outer, RSP_nz, RSP_T_anchor, RSP_T_inner, RSP_Teff, RSP_mass, RSP_L, RSP_X, RSP_Z, RSP_T_inner_tolerance, &
    RSP_relax_max_tries, RSP_hydro_only, &
   RSP_tau_surf_for_atm_grey_with_kap, RSP_fixed_Psurf, RSP_Avel, RSP_Arnd, RSP_relax_adjust_inner_mass_distribution, &
   use_other_RSP_linear_analysis, RSP_use_atm_grey_with_kap_for_Psurf, &
   RSP_fraction_1st_overtone,RSP_fraction_2nd_overtone, RSP_testing, RSP_use_Prad_for_Psurf, RSP_map_zone_interval, &
   RSP_write_map, RSP_map_filename, RSP_map_history_filename, RSP_map_first_period, RSP_map_last_period, &
   use_other_RSP_build_model, RSP_Psurf, RSP_work_period, RSP_work_filename, RSP_nmodes, RSP_surface_tau, &
   set_RSP_Psurf_to_multiple_of_initial_P1, use_RSP_new_start_scheme, RSP_do_check_omega, RSP_report_check_omega_changes, &
    RSP_relax_initial_model, RSP_trace_RSP_build_model, &
   RSP_GREKM_avg_abs_limit, RSP_GREKM_avg_abs_frac_new, RSP_kap_density_factor, RSP_map_columns_filename, &
   RSP_relax_alfap_before_alfat, RSP_max_outer_dm_tries, RSP_max_inner_scale_tries, RSP_T_anchor_tolerance, &
    ! mass gain or loss  
    mass_change, mass_change_full_on_dt, mass_change_full_off_dt, trace_dt_control_mass_change, &
    min_wind, max_wind, use_accreted_material_j, accreted_material_j, D_omega_mixing_rate, &
    D_omega_mixing_across_convection_boundary, max_q_for_D_omega_zero_in_convection_region, nu_omega_mixing_rate, &
    nu_omega_mixing_across_convection_boundary, max_q_for_nu_omega_zero_in_convection_region, &
    mdot_omega_power, max_rotational_mdot_boost, max_mdot_jump_for_rotation, &
    lim_trace_rotational_mdot_boost, rotational_mdot_boost_fac, rotational_mdot_kh_fac, surf_avg_tau, surf_avg_tau_min, &
    max_tries_for_implicit_wind, iwind_tolerance, iwind_lambda, super_eddington_scaling_factor, &
    super_eddington_wind_Ledd_factor, wind_boost_full_off_L_div_Ledd, wind_boost_full_on_L_div_Ledd, &
    super_eddington_wind_max_boost, trace_super_eddington_wind_boost, &
    rlo_scaling_factor, rlo_wind_min_L, rlo_wind_max_Teff, rlo_wind_roche_lobe_radius, &
    roche_lobe_xfer_full_on, roche_lobe_xfer_full_off, rlo_wind_base_mdot, rlo_wind_scale_height, &
    hot_wind_scheme, cool_wind_RGB_scheme, cool_wind_AGB_scheme, RGB_to_AGB_wind_switch, &
    Reimers_scaling_factor, Blocker_scaling_factor, de_Jager_scaling_factor, van_Loon_scaling_factor, &
    Nieuwenhuijzen_scaling_factor, Vink_scaling_factor, &
    Dutch_scaling_factor, Dutch_wind_lowT_scheme, wind_He_layer_limit, &
    wind_H_envelope_limit, wind_H_He_envelope_limit, hot_wind_full_on_T, cool_wind_full_on_T, &
    
    ! composition of added mass
    accrete_same_as_surface, &
    accrete_given_mass_fractions, num_accretion_species, accretion_species_id, accretion_species_xa, &
    accretion_h1, accretion_h2, accretion_he3, accretion_he4, accretion_zfracs, accretion_dump_missing_metals_into_heaviest, &
    
    ! special list of z fractions
    z_fraction_li, z_fraction_be, z_fraction_b, z_fraction_c, z_fraction_n,&
    z_fraction_o, z_fraction_f, z_fraction_ne, z_fraction_na, z_fraction_mg, z_fraction_al, &
    z_fraction_si, z_fraction_p, z_fraction_s, z_fraction_cl, z_fraction_ar, z_fraction_k, &
    z_fraction_ca, z_fraction_sc, z_fraction_ti, z_fraction_v, z_fraction_cr, z_fraction_mn, &
    z_fraction_fe, z_fraction_co, z_fraction_ni, z_fraction_cu, z_fraction_zn, &
    lgT_lo_for_set_new_abundances, lgT_hi_for_set_new_abundances, &
    
    ! automatic stops for mass loss/gain
    max_star_mass_for_gain, min_star_mass_for_loss, max_T_center_for_any_mass_loss, max_T_center_for_full_mass_loss, &
    
    ! extra power source
    extra_power_source, &
    
    ! relaxation parameters
    relax_dlnZ, relax_dY, &
    
    ! mesh adjustment
    show_mesh_changes, okay_to_remesh, restore_mesh_on_retry, num_steps_to_hold_mesh_after_retry, &
    max_rel_delta_IE_for_mesh_total_energy_balance, &
    trace_mesh_adjust_error_in_conservation, max_allowed_nz, mesh_max_allowed_ratio, &
    remesh_max_allowed_logT, max_delta_x_for_merge, &
    mesh_ok_to_merge, mesh_max_k_old_for_split, mesh_min_k_old_for_split, &
    mesh_adjust_get_T_from_E, &
    max_dq, min_dq, min_dq_for_split, min_dq_for_xa, min_dq_for_xa_convective, &
    mesh_min_dlnR, merge_if_dlnR_too_small, min_dq_for_logT, &
    mesh_min_dr_div_dRstar, merge_if_dr_div_dRstar_too_small, &
    mesh_min_dr_div_cs, merge_if_dr_div_cs_too_small, &
    max_center_cell_dq, max_surface_cell_dq, max_num_subcells, max_num_merge_cells, &
    mesh_delta_coeff, mesh_delta_coeff_for_highT, &
    logT_max_for_standard_mesh_delta_coeff, logT_min_for_highT_mesh_delta_coeff, remesh_dt_limit, &
    mesh_Pgas_div_P_exponent, &
    E_function_weight, E_function_param, P_function_weight, &
    mesh_logX_species, &
    mesh_logX_min_for_extra, mesh_dlogX_dlogP_extra, mesh_dlogX_dlogP_full_on, mesh_dlogX_dlogP_full_off, &

    mesh_dlog_eps_min_for_extra, mesh_dlog_eps_dlogP_full_on, mesh_dlog_eps_dlogP_full_off, &
    mesh_dlog_pp_dlogP_extra, mesh_dlog_cno_dlogP_extra, mesh_dlog_3alf_dlogP_extra, &
    mesh_dlog_burn_c_dlogP_extra, mesh_dlog_burn_n_dlogP_extra, mesh_dlog_burn_o_dlogP_extra, &
    mesh_dlog_burn_ne_dlogP_extra, mesh_dlog_burn_na_dlogP_extra, mesh_dlog_burn_mg_dlogP_extra, &
    mesh_dlog_burn_si_dlogP_extra, mesh_dlog_burn_s_dlogP_extra, mesh_dlog_burn_ar_dlogP_extra, &
    mesh_dlog_burn_ca_dlogP_extra, mesh_dlog_burn_ti_dlogP_extra, mesh_dlog_burn_cr_dlogP_extra, &
    mesh_dlog_burn_fe_dlogP_extra, mesh_dlog_cc_dlogP_extra, mesh_dlog_co_dlogP_extra, mesh_dlog_oo_dlogP_extra, &
    mesh_dlog_pnhe4_dlogP_extra, mesh_dlog_photo_dlogP_extra, mesh_dlog_other_dlogP_extra, &
    T_function1_weight, T_function2_weight, T_function2_param, mesh_delta_coeff_factor_smooth_iters, &
    R_function_weight, R_function_param, &
    R_function2_weight, R_function2_param1, R_function2_param2, &
    R_function3_weight, M_function_weight, M_function_param, &
    gradT_function_weight, log_tau_function_weight, log_kap_function_weight, omega_function_weight, &
    gam_function_weight, gam_function_param1, gam_function_param2, &
    xa_function_species, xa_function_weight, xa_function_param, xa_mesh_delta_coeff, split_merge_amr_mesh_delta_coeff, &
    use_split_merge_amr, split_merge_amr_nz_baseline, split_merge_amr_log_zoning, split_merge_amr_hybrid_zoning, &
    split_merge_amr_flipped_hybrid_zoning, split_merge_amr_logtau_zoning, split_merge_amr_okay_to_split_nz, split_merge_amr_nz_r_core, &
    split_merge_amr_okay_to_split_1, merge_amr_inhibit_at_jumps, split_merge_amr_MaxLong, split_merge_amr_nz_r_core_fraction, &
    split_merge_amr_MaxShort, merge_amr_max_abs_du_div_cs, &
    merge_amr_ignore_surface_cells, merge_amr_k_for_ignore_surface_cells, &
    merge_amr_du_div_cs_limit_only_for_compression, split_merge_amr_avoid_repeated_remesh, split_merge_amr_r_core_cm, &
    split_merge_amr_dq_min, split_merge_amr_dq_max, split_merge_amr_max_iters, trace_split_merge_amr, equal_split_density_amr, &

    ! nuclear reaction parameters
    screening_mode, default_net_name, net_logTcut_lo, net_logTcut_lim, &
    eps_nuc_factor, op_split_burn_eps_nuc_infall_limit, eps_WD_sedimentation_factor, &
    max_abs_eps_nuc, dxdt_nuc_factor, max_abar_for_burning, mix_factor, &
    fe56ec_fake_factor, min_T_for_fe56ec_fake_factor, weak_rate_factor, &
    sig_term_limit, sig_min_factor_for_high_Tcenter, Tcenter_min_for_sig_min_factor_full_on,&
    Tcenter_max_for_sig_min_factor_full_off, max_delta_m_to_bdy_for_sig_min_factor, &
    delta_m_lower_for_sig_min_factor, delta_m_upper_for_sig_min_factor, &
    am_sig_term_limit, am_D_mix_factor, am_gradmu_factor, am_nu_factor, &
    D_visc_factor, D_DSI_factor, D_SH_factor, D_SSI_factor, D_ES_factor, D_GSF_factor, D_ST_factor, &
    am_nu_non_rotation_factor, skip_rotation_in_convection_zones, am_nu_DSI_factor, am_nu_SH_factor,&
    am_nu_SSI_factor, am_nu_ES_factor, am_nu_GSF_factor, am_nu_ST_factor, am_nu_visc_factor, smooth_am_nu_rot, &
    am_nu_omega_rot_factor, am_nu_omega_non_rot_factor, am_nu_j_rot_factor, am_nu_j_non_rot_factor, &
    smooth_nu_ST, smooth_D_ST, smooth_D_SH, smooth_D_DSI, smooth_D_ES, smooth_D_SSI, smooth_D_GSF, smooth_D_omega, &
    do_adjust_J_lost, premix_omega, angular_momentum_error_warn, angular_momentum_error_retry, &
    simple_i_rot_flag, recalc_mixing_info_each_substep, adjust_J_fraction, &
    min_q_for_adjust_J_lost, min_J_div_delta_J, max_mdot_redo_cnt, mdot_revise_factor, &
    implicit_mdot_boost, min_years_dt_for_redo_mdot, surf_omega_div_omega_crit_limit, surf_omega_div_omega_crit_tol, &
    fitted_fp_ft_i_rot, w_div_wcrit_max, w_div_wcrit_max2, &
    fp_min, ft_min, fp_error_limit, ft_error_limit, &
    D_mix_rotation_max_logT_full_on, D_mix_rotation_min_logT_full_off, &
    set_uniform_am_nu_non_rot, uniform_am_nu_non_rot, &
    set_min_am_nu_non_rot, min_am_nu_non_rot, min_center_Ye_for_min_am_nu_non_rot, &
    set_min_D_mix_below_Tmax, min_D_mix_below_Tmax, set_min_D_mix_in_H_He, min_D_mix_in_H_He, &
    set_min_D_mix, mass_lower_limit_for_min_D_mix, mass_upper_limit_for_min_D_mix, &
    min_D_mix, min_center_Ye_for_min_D_mix, &
    smooth_outer_xa_big, smooth_outer_xa_small, nonlocal_NiCo_kap_gamma, nonlocal_NiCo_decay_heat, &
    dtau_gamma_NiCo_decay_heat, max_logT_for_net, reaction_neuQs_factor, &
    
    ! element diffusion parameters
    diffusion_use_iben_macdonald, diffusion_use_paquette, diffusion_use_cgs_solver, &
    diffusion_use_full_net, do_WD_sedimentation_heating, min_xa_for_WD_sedimentation_heating, &
    do_diffusion_heating, do_element_diffusion, &
    cgs_thermal_diffusion_eta_full_on, cgs_thermal_diffusion_eta_full_off, diffusion_min_dq_at_surface, &
    diffusion_min_T_at_surface, diffusion_min_dq_ratio_at_surface, diffusion_dt_limit, &
    diffusion_min_X_hard_limit, diffusion_X_total_atol, diffusion_X_total_rtol, &
    diffusion_upwind_abs_v_limit, diffusion_dt_div_timescale, diffusion_min_num_substeps, &
    diffusion_max_iters_per_substep, diffusion_max_retries_per_substep, diffusion_v_max, &
    diffusion_gamma_full_off, diffusion_gamma_full_on, diffusion_T_full_on, diffusion_T_full_off, &
    D_mix_ignore_diffusion, diffusion_calculates_ionization, diffusion_nsmooth_typical_charge, &
    diffusion_tol_correction_max, diffusion_tol_correction_norm, &
    diffusion_AD_dm_full_on, diffusion_AD_dm_full_off, diffusion_AD_boost_factor, &
    diffusion_SIG_factor, diffusion_GT_factor, &
    diffusion_Vlimit_dm_full_on, diffusion_Vlimit_dm_full_off, diffusion_Vlimit, &
    diffusion_max_T_for_radaccel, diffusion_min_T_for_radaccel, diffusion_max_Z_for_radaccel, &
    diffusion_min_Z_for_radaccel, diffusion_screening_for_radaccel, &
    op_mono_data_path, op_mono_data_cache_filename, &
    show_diffusion_info, show_diffusion_substep_info, show_diffusion_timing, &
    diffusion_num_classes, diffusion_class_representative, diffusion_class_A_max, &
    diffusion_class_typical_charge, diffusion_class_factor, &
    diffusion_use_isolve, diffusion_rtol_for_isolve, diffusion_atol_for_isolve, &
    diffusion_maxsteps_for_isolve, diffusion_isolve_solver, &
    
    ! eos controls
    use_fixed_XZ_for_eos, fixed_X_for_eos, fixed_Z_for_eos, use_d_eos_dxa, &
    report_eos_settings_at_start_of_run, &

    ! opacity controls
    use_simple_es_for_kap, use_starting_composition_for_kap, &
    min_kap_for_dPrad_dm_eqn, low_logT_op_mono_full_off, low_logT_op_mono_full_on, high_logT_op_mono_full_off, &
    high_logT_op_mono_full_on, op_mono_min_X_to_include, use_op_mono_alt_get_kap, &
    
    
    include_L_in_correction_limits, include_v_in_correction_limits, include_u_in_correction_limits, include_w_in_correction_limits, &
    
    ! asteroseismology controls
    get_delta_nu_from_scaled_solar, nu_max_sun, delta_nu_sun, Teff_sun, delta_Pg_mode_freq, &
    
    ! hydro parameters
    opacity_factor, opacity_max, min_logT_for_opacity_factor_off, min_logT_for_opacity_factor_on, &
    max_logT_for_opacity_factor_on, max_logT_for_opacity_factor_off, &
    non_nuc_neu_factor, always_use_dedt_form_of_energy_eqn, &
    use_dedt_form_of_energy_eqn, use_time_centered_eps_grav, &
    use_mass_corrections, use_gravity_rotation_correction, eps_grav_factor, eps_mdot_factor, &
    include_composition_in_eps_grav, no_dedt_form_during_relax, &
    max_abs_rel_change_surf_lnS, always_use_eps_grav_form_of_energy_eqn, &
    max_num_surf_revisions, Gamma_lnS_eps_grav_full_off, Gamma_lnS_eps_grav_full_on, &
    use_dPrad_dm_form_of_T_gradient_eqn, use_gradT_actual_vs_gradT_MLT_for_T_gradient_eqn, dedt_eqn_r_scale, &
    RTI_A, RTI_B, RTI_C, RTI_D, RTI_max_alpha, RTI_C_X_factor, RTI_C_X0_frac, steps_before_use_velocity_time_centering, &
    RTI_dm_for_center_eta_nondecreasing, RTI_min_dm_behind_shock_for_full_on, RTI_energy_floor, &
    RTI_D_mix_floor, RTI_min_m_for_D_mix_floor, RTI_log_max_boost, RTI_m_full_boost, RTI_m_no_boost, &
    conv_vel_D, conv_vel_siglimit, conv_vel_v0, include_P_in_velocity_time_centering, include_L_in_velocity_time_centering, &
    P_theta_for_velocity_time_centering, L_theta_for_velocity_time_centering, &
    min_q_for_normal_mlt_gradT_full_off, max_q_for_normal_mlt_gradT_full_on, &
    conv_vel_ignore_thermohaline, conv_vel_ignore_semiconvection, use_P_d_1_div_rho_form_of_work_when_time_centering_velocity, &
    conv_vel_fully_lagrangian, conv_vel_include_homologous_term, conv_vel_use_mlt_vc_start, compare_TDC_to_MLT, &
    velocity_logT_lower_bound, max_dt_yrs_for_velocity_logT_lower_bound, velocity_q_upper_bound, &

    ! hydro solver
    use_gold2_tolerances, gold2_solver_iters_timestep_limit, steps_before_use_gold2_tolerances, &
    gold2_iter_for_resid_tol2, gold2_iter_for_resid_tol3, gold2_tol_residual_norm1, gold2_tol_max_residual1, &
    gold2_tol_residual_norm2, gold2_tol_max_residual2, gold2_tol_residual_norm3, gold2_tol_max_residual3, &
    tol_correction_norm, tol_max_correction, correction_xa_limit, &
    tol_correction_high_T_limit, tol_correction_norm_high_T, tol_max_correction_high_T, &
    tol_correction_extreme_T_limit, tol_correction_norm_extreme_T, tol_max_correction_extreme_T, &
    tol_bad_max_correction, bad_max_correction_series_limit, &
    tol_max_residual1, tol_residual_norm1, tol_max_residual2, &
    tol_residual_norm2, tol_max_residual3, tol_residual_norm3, &
    warning_limit_for_max_residual, trace_solver_damping, &
    relax_use_gold_tolerances, relax_tol_correction_norm, relax_tol_max_correction, &
    relax_solver_iters_timestep_limit, &
    relax_iter_for_resid_tol2, relax_tol_residual_norm1, relax_tol_max_residual1, &
    relax_iter_for_resid_tol3, relax_tol_residual_norm2, relax_tol_max_residual2, &
    relax_tol_residual_norm3, relax_tol_max_residual3, relax_maxT_for_gold_tolerances, &
    use_gold_tolerances, gold_solver_iters_timestep_limit, maxT_for_gold_tolerances, &
    gold_tol_residual_norm1, gold_tol_max_residual1, gold_iter_for_resid_tol2, &
    gold_tol_residual_norm2, gold_tol_max_residual2, gold_iter_for_resid_tol3, &
    gold_tol_residual_norm3, gold_tol_max_residual3 , steps_before_use_gold_tolerances, &
    include_rotation_in_total_energy, convergence_ignore_equL_residuals, convergence_ignore_alpha_RTI_residuals, &
    iter_for_resid_tol2, iter_for_resid_tol3, &
    solver_itermin, solver_itermin_until_reduce_min_corr_coeff, &
    solver_reduced_min_corr_coeff, do_solver_damping_for_neg_xa, &
    hydro_mtx_max_allowed_abs_dlogT, hydro_mtx_max_allowed_abs_dlogRho, &
    min_logT_for_hydro_mtx_max_allowed, hydro_mtx_max_allowed_logT, &
    hydro_mtx_max_allowed_logRho, report_min_rcond_from_DGESXV, &
    hydro_mtx_min_allowed_logT, hydro_mtx_min_allowed_logRho, use_DGESVX_in_bcyclic, use_equilibration_in_DGESVX, &
    op_split_burn, op_split_burn_min_T, op_split_burn_eps, op_split_burn_odescal, &
    op_split_burn_min_T_for_variable_T_solver, solver_test_partials_show_dx_var_name, &
    tiny_corr_coeff_limit, scale_correction_norm, corr_param_factor, num_times_solver_reuse_mtx, &
    scale_max_correction, ignore_min_corr_coeff_for_scale_max_correction, ignore_too_large_correction, ignore_species_in_max_correction, &
    corr_norm_jump_limit, max_corr_jump_limit, resid_norm_jump_limit, max_resid_jump_limit, RSP2_use_mass_interp_face_values, &
    corr_coeff_limit, tiny_corr_factor, solver_test_partials_call_number, solver_test_partials_iter_number, &
    max_tries1, solver_max_tries_before_reject, max_tries_for_retry, max_tries_after_5_retries, solver_test_partials_sink_name, &
    max_tries_after_10_retries, max_tries_after_20_retries, retry_limit, redo_limit, use_Pvsc_art_visc, Pvsc_cq, Pvsc_zsh, &
    min_xa_hard_limit, min_xa_hard_limit_for_highT, logT_max_for_min_xa_hard_limit, logT_min_for_min_xa_hard_limit_for_highT, &
    sum_xa_hard_limit, sum_xa_hard_limit_for_highT, logT_max_for_sum_xa_hard_limit, logT_min_for_sum_xa_hard_limit_for_highT, &
    xa_clip_limit, report_solver_progress, solver_test_partials_k_high, RSP2_use_L_eqn_at_surface, RSP2_use_RSP_eqn_for_Y_face, &
    solver_epsder_chem, solver_epsder_struct, solver_numerical_jacobian, energy_conservation_dump_model_number, &
    solver_jacobian_nzlo, solver_jacobian_nzhi, solver_check_everything, solver_inspect_soln_flag, RSP2_assume_HSE, &
    solver_test_partials_dx_0, solver_test_partials_k, solver_show_correction_info, eps_mdot_leak_frac_factor, &
    solver_test_partials_write_eos_call_info, solver_save_photo_call_number, RSP2_min_Lc_div_L_for_convective_mixing_type, &
    solver_test_partials_var_name, solver_test_partials_equ_name, RSP2_min_Lt_div_L_for_overshooting_mixing_type, &
    solver_test_eos_partials, solver_test_kap_partials, solver_test_net_partials, solver_test_atm_partials, &
    fill_arrays_with_NaNs, zero_when_allocate, warn_when_large_rel_run_E_err, solver_test_partials_k_low, &
    warn_when_large_virial_thm_rel_err, warn_when_get_a_bad_eos_result, warn_rates_for_high_temp, max_safe_logT_for_rates, &
    RSP2_alfap, RSP2_alfat, RSP2_alfam, RSP2_alfar, RSP2_Lsurf_factor, RSP2_use_Stellingwerf_Lr, &
    RSP2_alfad, RSP2_num_outermost_cells_forced_nonturbulent, RSP2_num_innermost_cells_forced_nonturbulent, &
    RSP2_target_steps_per_cycle, RSP2_max_num_periods, RSP2_work_period, RSP2_map_first_period, RSP2_map_last_period, &
    RSP2_min_max_R_for_periods, RSP2_GREKM_avg_abs_frac_new, RSP2_GREKM_avg_abs_limit, RSP2_map_zone_interval, &
    RSP2_work_filename, RSP2_map_columns_filename, RSP2_map_filename, RSP2_map_history_filename, RSP2_write_map, &
    RSP2_min_dt_div_tau_conv_switch_to_MLT, RSP2_min_dt_years_switch_to_MLT, &
    RSP2_w_min_for_damping, RSP2_source_seed, RSP2_w_fix_if_neg, max_X_for_conv_timescale, min_X_for_conv_timescale, &
    max_q_for_conv_timescale, min_q_for_conv_timescale, max_q_for_QHSE_timescale, min_q_for_QHSE_timescale, &
    
    
    ! timestep
    time_delta_coeff, min_timestep_factor, max_timestep_factor, timestep_factor_for_retries, retry_hold, &
    neg_mass_fraction_hold, timestep_dt_factor, use_dt_low_pass_controller, &
    force_timestep_min, force_timestep_min_years, force_timestep_min_factor, force_timestep, force_timestep_years, &
    varcontrol_target, min_allowed_varcontrol_target, varcontrol_dt_limit_ratio_hard_max, xa_scale, &
    solver_iters_timestep_limit, burn_steps_limit, burn_steps_hard_limit, &
    diffusion_steps_limit, diffusion_steps_hard_limit, diffusion_iters_limit, diffusion_iters_hard_limit, &
    dt_div_dt_cell_collapse_limit, dt_div_dt_cell_collapse_hard_limit, &
    dt_div_min_dr_div_cs_limit, dt_div_min_dr_div_cs_hard_limit, min_abs_u_div_cs_for_dt_div_min_dr_div_cs_limit, &
    min_abs_du_div_cs_for_dt_div_min_dr_div_cs_limit, min_k_for_dt_div_min_dr_div_cs_limit, &
    min_q_for_dt_div_min_dr_div_cs_limit, max_q_for_dt_div_min_dr_div_cs_limit, check_remnant_only_for_dt_div_min_dr_div_cs_limit, &
    dX_mix_dist_limit, dH_limit_min_H, dH_limit, dH_hard_limit, dH_div_H_limit_min_H, &
    dH_div_H_limit, dH_div_H_hard_limit, dH_decreases_only, max_timestep_factor_at_high_T, &
    dHe_limit_min_He, dHe_limit, dHe_hard_limit, dHe_div_He_limit_min_He, &
    dHe_div_He_limit, dHe_div_He_hard_limit, dHe_decreases_only, min_logT_for_max_timestep_factor_at_high_T, &
    dHe3_limit_min_He3, dHe3_limit, dHe3_hard_limit, dHe3_div_He3_limit_min_He3, &
    dHe3_div_He3_limit, dHe3_div_He3_hard_limit, dHe3_decreases_only, dX_div_X_at_high_T_limit_lgT_min, &
    dX_limit_min_X, dX_limit, dX_hard_limit, dX_div_X_limit_min_X, dX_div_X_at_high_T_hard_limit, &
    dX_div_X_limit, dX_div_X_hard_limit, dX_decreases_only, dX_div_X_at_high_T_limit, &
    dX_nuc_drop_min_X_limit, dX_nuc_drop_max_A_limit, dX_nuc_drop_limit, &
    dX_nuc_drop_limit_at_high_T, dX_nuc_drop_hard_limit, dX_nuc_drop_min_yrs_for_dt, &
    dL_div_L_limit_min_L, dL_div_L_limit, dL_div_L_hard_limit, &
    delta_lgP_limit, delta_lgP_hard_limit, delta_lgP_limit_min_lgP, &
    delta_lgRho_limit, delta_lgRho_hard_limit, delta_lgRho_limit_min_lgRho, &
    delta_lgT_limit, delta_lgT_hard_limit, delta_lgT_limit_min_lgT, &
    delta_lgE_limit, delta_lgE_hard_limit, delta_lgE_limit_min_lgE, &
    delta_lgR_limit, delta_lgR_hard_limit, delta_lgR_limit_min_lgR, &
    delta_Ye_highT_limit, delta_Ye_highT_hard_limit, minT_for_highT_Ye_limit, &
    delta_lgL_nuc_cat_limit, delta_lgL_nuc_cat_hard_limit, lgL_nuc_cat_burn_min, lgL_nuc_mix_dist_limit, &
    delta_lgT_max_limit_only_after_near_zams, &
    delta_lgL_H_limit, delta_lgL_H_hard_limit, lgL_H_burn_min, lgL_H_drop_factor, lgL_H_burn_relative_limit, &
    delta_lgL_He_limit, delta_lgL_He_hard_limit, lgL_He_burn_min, lgL_He_drop_factor, lgL_He_burn_relative_limit, &
    delta_lgL_z_limit, delta_lgL_z_hard_limit, lgL_z_burn_min, lgL_z_drop_factor, lgL_z_burn_relative_limit, &
    delta_lgL_power_photo_limit, delta_lgL_power_photo_hard_limit, lgL_power_photo_burn_min, lgL_power_photo_drop_factor, &
    delta_lgL_nuc_limit, delta_lgL_nuc_hard_limit, lgL_nuc_burn_min, lgL_nuc_drop_factor, min_lgT_for_lgL_power_photo_limit, &
    delta_lgL_nuc_at_high_T_limit, delta_lgL_nuc_at_high_T_hard_limit, delta_lgL_nuc_at_high_T_limit_lgT_min, &
    delta_lgRho_cntr_limit, delta_lgRho_cntr_hard_limit, delta_lgP_cntr_limit, delta_lgP_cntr_hard_limit, &
    delta_lgT_cntr_limit, delta_lgT_cntr_hard_limit, delta_lgT_cntr_limit_only_after_near_zams, &
    delta_lgT_max_limit, delta_lgT_max_hard_limit, max_lgT_for_lgL_nuc_limit, &
    delta_log_eps_nuc_limit, delta_log_eps_nuc_hard_limit, delta_lgT_max_limit_lgT_min, &
    delta_lgT_max_at_high_T_limit, delta_lgT_max_at_high_T_hard_limit, delta_lgT_max_at_high_T_limit_lgT_min, &
    delta_dX_div_X_cntr_min, delta_dX_div_X_cntr_max, delta_dX_div_X_cntr_limit, delta_dX_div_X_cntr_hard_limit, &
    delta_dX_div_X_drop_only, delta_lg_XH_drop_only, &
    delta_lg_XHe_drop_only, delta_lg_XC_drop_only, delta_lg_XNe_drop_only, delta_lg_XO_drop_only, delta_lg_XSi_drop_only, &
    delta_XH_drop_only, delta_XHe_drop_only, delta_XC_drop_only, delta_XNe_drop_only, delta_XO_drop_only, delta_XSi_drop_only, &    
    delta_lg_XH_cntr_min, delta_lg_XH_cntr_max, delta_lg_XH_cntr_limit, delta_lg_XH_cntr_hard_limit, &
    delta_lg_XHe_cntr_min, delta_lg_XHe_cntr_max, delta_lg_XHe_cntr_limit, delta_lg_XHe_cntr_hard_limit, &
    delta_lg_XC_cntr_min, delta_lg_XC_cntr_max, delta_lg_XC_cntr_limit, delta_lg_XC_cntr_hard_limit, &
    delta_lg_XNe_cntr_limit, delta_lg_XNe_cntr_hard_limit, delta_lg_XNe_cntr_min, delta_lg_XNe_cntr_max, &
    delta_lg_XO_cntr_limit, delta_lg_XO_cntr_hard_limit, delta_lg_XO_cntr_min, delta_lg_XO_cntr_max, &
    delta_lg_XSi_cntr_limit, delta_lg_XSi_cntr_hard_limit, delta_lg_XSi_cntr_min, delta_lg_XSi_cntr_max, &
    delta_XH_cntr_limit, delta_XH_cntr_hard_limit, delta_XHe_cntr_limit, delta_XHe_cntr_hard_limit, &
    delta_XC_cntr_limit, delta_XC_cntr_hard_limit, delta_XNe_cntr_limit, delta_XNe_cntr_hard_limit, &
    delta_XO_cntr_limit, delta_XO_cntr_hard_limit, delta_XSi_cntr_limit, delta_XSi_cntr_hard_limit, &
    delta_lgTeff_limit, delta_lgTeff_hard_limit, delta_lgL_limit, delta_lgL_limit_L_min, delta_lgL_hard_limit, &
    delta_HR_ds_L, delta_HR_ds_Teff, delta_HR_limit, delta_HR_hard_limit, &
    delta_lg_star_mass_limit, delta_lg_star_mass_hard_limit, &
    delta_mdot_atol, delta_mdot_rtol, delta_mdot_limit, delta_mdot_hard_limit, &
    adjust_J_q_limit, adjust_J_q_hard_limit, &
    never_skip_hard_limits, relax_hard_limits_after_retry, &
    report_dt_hard_limit_retries, report_min_dr_div_cs, report_solver_dt_info, &
    limit_for_rel_error_in_energy_conservation, hard_limit_for_rel_error_in_energy_conservation, &

    ! atmosphere -- surface boundary conditions

    atm_option, atm_off_table_option, Pextra_factor, &
    atm_fixed_Teff, atm_fixed_Psurf, atm_fixed_Tsurf, &

    atm_T_tau_relation, atm_T_tau_opacity, atm_T_tau_errtol, atm_T_tau_max_iters, &
    atm_T_tau_max_steps, &

    atm_table, &

    atm_irradiated_opacity, atm_irradiated_errtol, atm_irradiated_T_eq, &
    atm_irradiated_kap_v, atm_irradiated_kap_v_div_kap_th, atm_irradiated_P_surf, &
    atm_irradiated_max_iters, &

    use_compression_outer_BC, use_momentum_outer_BC, Tsurf_factor, use_zero_Pgas_outer_BC, &
    fixed_Psurf, use_fixed_Psurf_outer_BC, fixed_vsurf, use_fixed_vsurf_outer_BC, &
    
    atm_build_tau_outer, atm_build_dlogtau, atm_build_errtol, &

    use_T_tau_gradr_factor, &
    
    ! extra heat near surface to model irradiation
    irradiation_flux, column_depth_for_irradiation, &
    
    ! uniform extra heat
    inject_uniform_extra_heat, min_q_for_uniform_extra_heat, max_q_for_uniform_extra_heat, &
    inject_extra_ergs_sec, base_of_inject_extra_ergs_sec, total_mass_for_inject_extra_ergs_sec, &
    start_time_for_inject_extra_ergs_sec, duration_for_inject_extra_ergs_sec, &
    inject_until_reach_model_with_total_energy, &
    
    ! mass gain or loss
    no_wind_if_no_rotation, max_logT_for_k_below_const_q, &
    max_q_for_k_below_const_q, min_q_for_k_below_const_q, max_logT_for_k_const_mass, &
    min_q_for_k_const_mass, max_q_for_k_const_mass, &
    
    ! info for debugging
    stop_for_bad_nums, report_ierr, report_bad_negative_xa, diffusion_dump_call_number, &
     
    ! controls for the evolve routine
    trace_evolve, &

    ! misc
    min_chem_eqn_scale, zams_filename, set_rho_to_dm_div_dV, use_other_momentum_implicit, &
    use_other_eos, use_other_surface_PT, use_other_kap, use_other_diffusion, use_other_diffusion_factor, &
    use_other_adjust_mdot, use_other_j_for_adjust_J_lost, use_other_alpha_mlt, use_other_remove_surface, &
    use_other_am_mixing, use_other_brunt, use_other_brunt_smoothing, use_other_solver_monitor, &
    use_other_build_initial_model, use_other_cgrav, use_other_energy_implicit, use_other_momentum, &
    use_other_energy, use_other_mesh_functions, use_other_eps_grav, use_other_gradr_factor, &
    use_other_D_mix, use_other_neu, use_other_net_get, use_other_opacity_factor, use_other_pressure, &
    use_other_diffusion_coefficients, use_other_pgstar_plots, use_other_eval_fp_ft, use_other_eval_i_rot, use_other_torque, &
    use_other_torque_implicit, use_other_wind, use_other_accreting_state, use_other_after_struct_burn_mix, use_other_mesh_delta_coeff_factor, &
    use_other_before_struct_burn_mix, use_other_astero_freq_corr, use_other_timestep_limit, use_other_set_pgstar_controls, &
    use_other_screening, &
    x_ctrl, x_integer_ctrl, x_logical_ctrl, x_character_ctrl, &
    
    ! extra files
    read_extra_controls_inlist1, extra_controls_inlist1_name, read_extra_controls_inlist2, &
    extra_controls_inlist2_name, read_extra_controls_inlist3, extra_controls_inlist3_name, &
    read_extra_controls_inlist4, extra_controls_inlist4_name, read_extra_controls_inlist5, extra_controls_inlist5_name, &
    save_controls_namelist, controls_namelist_name


 contains


 subroutine write_controls(s, fname, ierr)
    type (star_info), pointer :: s
 character (len=*), intent(in) :: fname
 integer, intent(out) :: ierr 
 integer :: iounit
 character (len=256) :: filename
 
 ierr = 0
 filename = fname

 if (len_trim(filename) == 0) filename = 'dump_controls.txt'

 open(newunit=iounit, file=trim(filename), &
    action='write', status='replace', iostat=ierr)
 if (ierr /= 0) then
    write(*,*) 'failed to open ' // trim(filename)
    return
 endif
 
 call set_controls_for_writing(s, ierr)
 if (ierr /= 0) then
    close(iounit)
    return
 end if

 write(iounit, nml=controls, iostat=ierr)
 
 write(*,*) 'write controls namelist values to "' // trim(filename)//'"'
    
 close(iounit)

 end subroutine write_controls


 subroutine do_one_setup(id, inlist, ierr)
    use utils_lib
    character (len=*), intent(in) :: inlist
    integer, intent(in) :: id
    integer, intent(out) :: ierr
    type (star_info), pointer :: s
    include 'formats'
    call get_star_ptr(id, s, ierr)
    if (ierr /= 0) return
    call set_default_controls
    call read_controls(id, inlist, ierr)
    if (ierr /= 0) return

   if (save_controls_namelist) then
      call write_controls(s, controls_namelist_name, ierr)
      if (ierr /= 0) then
         write(*,*) 'ierr from write_controls ' // &
            trim(controls_namelist_name)
         return
      end if
   end if

 end subroutine do_one_setup


 subroutine read_controls(id, filename, ierr)
 use star_private_def
 use utils_lib
 character(*), intent(in) :: filename
 integer, intent(in) :: id
 integer, intent(out) :: ierr

 type (star_info), pointer :: s
 ierr = 0
 call get_star_ptr(id, s, ierr)
 if (ierr /= 0) return

 call read_controls_file(s, filename, 1, ierr)
 call mkdir(s% photo_directory)
 call mkdir(s% log_directory)

 end subroutine read_controls


 recursive subroutine read_controls_file(s, filename, level, ierr)
 use star_private_def
 use utils_lib
 character(*), intent(in) :: filename
 type (star_info), pointer :: s
 integer, intent(in) :: level
 integer, intent(out) :: ierr
 logical :: read_extra1, read_extra2, read_extra3, read_extra4, read_extra5
 character (len=strlen) :: message, extra1, extra2, extra3, extra4, extra5
 integer :: unit

 ierr = 0

 if (level >= 10) then
 write(*,*) 'ERROR: too many levels of nested extra controls inlist files'
 ierr = -1
 return
 end if

 if (len_trim(filename) > 0) then
    open(newunit=unit, file=trim(filename), action='read', delim='quote', status='old', iostat=ierr)
    if (ierr /= 0) then
       write(*, *) 'Failed to open control namelist file ', trim(filename)
       return
    end if
    read(unit, nml=controls, iostat=ierr)
    close(unit)
    if (ierr /= 0) then
       write(*, *)
       write(*, *)
       write(*, *)
       write(*, *)
       write(*, '(a)') 'Failed while trying to read control namelist file: ' // trim(filename)
       write(*, '(a)') 'Perhaps the following runtime error message will help you find the problem.'
       write(*, *)
       open(newunit=unit, file=trim(filename), action='read', delim='quote', status='old', iostat=ierr)
       read(unit, nml=controls)
       close(unit)
       return
    end if
 end if

 call store_controls(s, ierr)

 ! recursive calls to read other inlists

 read_extra1 = read_extra_controls_inlist1
 read_extra_controls_inlist1 = .false.
 extra1 = extra_controls_inlist1_name
 extra_controls_inlist1_name = 'undefined'

 read_extra2 = read_extra_controls_inlist2
 read_extra_controls_inlist2 = .false.
 extra2 = extra_controls_inlist2_name
 extra_controls_inlist2_name = 'undefined'

 read_extra3 = read_extra_controls_inlist3
 read_extra_controls_inlist3 = .false.
 extra3 = extra_controls_inlist3_name
 extra_controls_inlist3_name = 'undefined'

 read_extra4 = read_extra_controls_inlist4
 read_extra_controls_inlist4 = .false.
 extra4 = extra_controls_inlist4_name
 extra_controls_inlist4_name = 'undefined'

 read_extra5 = read_extra_controls_inlist5
 read_extra_controls_inlist5 = .false.
 extra5 = extra_controls_inlist5_name
 extra_controls_inlist5_name = 'undefined'

 if (read_extra1) then
 write(*,*) 'read ' // trim(extra1)
 call read_controls_file(s, extra1, level+1, ierr)
 if (ierr /= 0) return
 end if

 if (read_extra2) then
 write(*,*) 'read ' // trim(extra2)
 call read_controls_file(s, extra2, level+1, ierr)
 if (ierr /= 0) return
 end if

 if (read_extra3) then
 write(*,*) 'read ' // trim(extra3)
 call read_controls_file(s, extra3, level+1, ierr)
 if (ierr /= 0) return
 end if

 if (read_extra4) then
    write(*,*) 'read ' // trim(extra4)
    call read_controls_file(s, extra4, level+1, ierr)
    if (ierr /= 0) return
 end if

 if (read_extra5) then
    write(*,*) 'read ' // trim(extra5)
    call read_controls_file(s, extra5, level+1, ierr)
    if (ierr /= 0) return
 end if

 end subroutine read_controls_file
 

 subroutine set_default_controls

 xa_central_lower_limit_species(:) = ''
 xa_central_lower_limit(:) = 0
 xa_central_upper_limit_species(:) = ''
 xa_central_upper_limit(:) = 0
 xa_surface_lower_limit_species(:) = ''
 xa_surface_lower_limit(:) = 0
 xa_surface_upper_limit_species(:) = ''
 xa_surface_upper_limit(:) = 0
 xa_average_lower_limit_species(:) = ''
 xa_average_lower_limit(:) = 0
 xa_average_upper_limit_species(:) = ''
 xa_average_upper_limit(:) = 0

 trace_history_value_name(:) = ''
 accretion_species_id(:) = ''
 accretion_species_xa(:) = 0
 mesh_logX_species(:) = ''
 mesh_logX_min_for_extra(:) = 0
 mesh_dlogX_dlogP_extra(:) = 0
 mesh_dlogX_dlogP_full_on(:) = 0
 mesh_dlogX_dlogP_full_off(:) = 0

 xa_function_species(:) = ''
 xa_function_weight(:) = 0
 xa_function_param(:) = 0
 xa_mesh_delta_coeff(:) = 0

 diffusion_class_representative(:) = ''
 diffusion_class_A_max(:) = 0
 diffusion_class_typical_charge(:) = 0
 diffusion_class_factor(:) = 0

 x_character_ctrl(:) = ''

 predictive_mix(:) = .FALSE.
 predictive_zone_type(:) = ''
 predictive_zone_loc(:) = ''
 predictive_bdy_loc(:) = ''
 predictive_bdy_q_min(:) = 0d0
 predictive_bdy_q_max(:) = 1d0
 predictive_superad_thresh(:) = 0d0
 predictive_avoid_reversal(:) = ''
 predictive_limit_ingestion(:) = ''
 predictive_ingestion_factor(:) = 0d0

 include 'controls.defaults'
 include 'controls_dev.defaults'

 end subroutine set_default_controls


 subroutine store_controls(s, ierr)
 use star_private_def
 use chem_def ! categories
 use utils_lib, only: mkdir
 type (star_info), pointer :: s
 integer, intent(out) :: ierr

 ierr = 0

 ! where to start
 s% initial_mass = initial_mass
 s% initial_z = initial_z
 s% initial_y = initial_y
 s% initial_he3 = initial_he3

 ! definition of core boundaries
 s% he_core_boundary_h1_fraction = he_core_boundary_h1_fraction
 s% co_core_boundary_he4_fraction = co_core_boundary_he4_fraction
 s% one_core_boundary_he4_c12_fraction = one_core_boundary_he4_c12_fraction
 s% fe_core_boundary_si28_fraction = fe_core_boundary_si28_fraction
 s% neutron_rich_core_boundary_Ye_max = neutron_rich_core_boundary_Ye_max
 s% min_boundary_fraction = min_boundary_fraction

 ! when to stop
 s% max_model_number = max_model_number
 s% max_number_retries = max_number_retries
 s% max_abs_rel_run_E_err = max_abs_rel_run_E_err
 s% relax_max_number_retries = relax_max_number_retries
 s% max_age = max_age
 s% max_age_in_days = max_age_in_days
 s% max_age_in_seconds = max_age_in_seconds
 s% num_adjusted_dt_steps_before_max_age = num_adjusted_dt_steps_before_max_age
 s% dt_years_for_steps_before_max_age = dt_years_for_steps_before_max_age
 s% reduction_factor_for_max_timestep = reduction_factor_for_max_timestep
 s% when_to_stop_rtol = when_to_stop_rtol
 s% when_to_stop_atol = when_to_stop_atol
 s% gamma_center_limit = gamma_center_limit
 s% eta_center_limit = eta_center_limit
 s% log_center_temp_limit = log_center_temp_limit
 s% log_max_temp_upper_limit = log_max_temp_upper_limit
 s% log_max_temp_lower_limit = log_max_temp_lower_limit
 s% log_center_temp_lower_limit = log_center_temp_lower_limit
 s% log_center_density_limit = log_center_density_limit
 s% log_center_density_lower_limit = log_center_density_lower_limit
 s% min_timestep_limit = min_timestep_limit

 s% center_entropy_limit = center_entropy_limit
 s% center_entropy_lower_limit = center_entropy_lower_limit
 s% max_entropy_limit = max_entropy_limit
 s% max_entropy_lower_limit = max_entropy_lower_limit

 s% fe_core_infall_limit = fe_core_infall_limit
 s% center_Ye_lower_limit = center_Ye_lower_limit
 s% center_R_lower_limit = center_R_lower_limit
 s% non_fe_core_infall_limit = non_fe_core_infall_limit
 s% non_fe_core_rebound_limit = non_fe_core_rebound_limit
 s% v_div_csound_surf_limit = v_div_csound_surf_limit
 s% v_div_csound_max_limit = v_div_csound_max_limit
 s% Lnuc_div_L_upper_limit = Lnuc_div_L_upper_limit
 s% Lnuc_div_L_lower_limit = Lnuc_div_L_lower_limit
 s% v_surf_div_v_kh_upper_limit = v_surf_div_v_kh_upper_limit
 s% v_surf_div_v_kh_lower_limit = v_surf_div_v_kh_lower_limit
 s% v_surf_div_v_esc_limit = v_surf_div_v_esc_limit
 s% v_surf_kms_limit = v_surf_kms_limit

 s% stop_near_zams = stop_near_zams
 s% stop_at_phase_PreMS = stop_at_phase_PreMS
 s% stop_at_phase_ZAMS = stop_at_phase_ZAMS
 s% stop_at_phase_IAMS = stop_at_phase_IAMS
 s% stop_at_phase_TAMS = stop_at_phase_TAMS
 s% stop_at_phase_He_Burn = stop_at_phase_He_Burn
 s% stop_at_phase_ZACHeB = stop_at_phase_ZACHeB
 s% stop_at_phase_TACHeB = stop_at_phase_TACHeB
 s% stop_at_phase_TP_AGB = stop_at_phase_TP_AGB
 s% stop_at_phase_C_Burn = stop_at_phase_C_Burn
 s% stop_at_phase_Ne_Burn = stop_at_phase_Ne_Burn
 s% stop_at_phase_O_Burn = stop_at_phase_O_Burn
 s% stop_at_phase_Si_Burn = stop_at_phase_Si_Burn
 s% stop_at_phase_WDCS = stop_at_phase_WDCS
 s% Lnuc_div_L_zams_limit = Lnuc_div_L_zams_limit
 s% gamma1_limit = gamma1_limit
 s% gamma1_limit_max_q = gamma1_limit_max_q
 s% gamma1_limit_max_v_div_vesc = gamma1_limit_max_v_div_vesc
 s% Pgas_div_P_limit = Pgas_div_P_limit
 s% Pgas_div_P_limit_max_q = Pgas_div_P_limit_max_q
 s% peak_burn_vconv_div_cs_limit = peak_burn_vconv_div_cs_limit
 s% omega_div_omega_crit_limit = omega_div_omega_crit_limit
 s% delta_nu_lower_limit = delta_nu_lower_limit
 s% delta_nu_upper_limit = delta_nu_upper_limit
 s% delta_Pg_lower_limit = delta_Pg_lower_limit
 s% delta_Pg_upper_limit = delta_Pg_upper_limit
 s% shock_mass_upper_limit = shock_mass_upper_limit
 s% mach1_mass_upper_limit = mach1_mass_upper_limit
 s% stop_when_reach_this_cumulative_extra_heating = stop_when_reach_this_cumulative_extra_heating

 s% xa_central_lower_limit_species = xa_central_lower_limit_species
 s% xa_central_lower_limit = xa_central_lower_limit

 s% xa_central_upper_limit_species = xa_central_upper_limit_species
 s% xa_central_upper_limit = xa_central_upper_limit

 s% xa_surface_lower_limit_species = xa_surface_lower_limit_species
 s% xa_surface_lower_limit = xa_surface_lower_limit

 s% xa_surface_upper_limit_species = xa_surface_upper_limit_species
 s% xa_surface_upper_limit = xa_surface_upper_limit

 s% xa_average_lower_limit_species = xa_average_lower_limit_species
 s% xa_average_lower_limit = xa_average_lower_limit

 s% xa_average_upper_limit_species = xa_average_upper_limit_species
 s% xa_average_upper_limit = xa_average_upper_limit

 s% HB_limit = HB_limit

 s% star_mass_max_limit = star_mass_max_limit
 s% star_mass_min_limit = star_mass_min_limit
 s% ejecta_mass_max_limit = ejecta_mass_max_limit
 s% remnant_mass_min_limit = remnant_mass_min_limit
 
 s% star_species_mass_min_limit = star_species_mass_min_limit
 s% star_species_mass_min_limit_iso = star_species_mass_min_limit_iso
 s% star_species_mass_max_limit = star_species_mass_max_limit
 s% star_species_mass_max_limit_iso = star_species_mass_max_limit_iso
 
 s% xmstar_min_limit = xmstar_min_limit
 s% xmstar_max_limit = xmstar_max_limit
 s% envelope_mass_limit = envelope_mass_limit
 s% envelope_fraction_left_limit = envelope_fraction_left_limit

 s% he_core_mass_limit = he_core_mass_limit
 s% co_core_mass_limit = co_core_mass_limit
 s% one_core_mass_limit = one_core_mass_limit
 s% fe_core_mass_limit = fe_core_mass_limit
 s% neutron_rich_core_mass_limit = neutron_rich_core_mass_limit

 s% he_layer_mass_lower_limit = he_layer_mass_lower_limit
 s% abs_diff_lg_LH_lg_Ls_limit = abs_diff_lg_LH_lg_Ls_limit
 s% Teff_upper_limit = Teff_upper_limit
 s% Teff_lower_limit = Teff_lower_limit
 s% photosphere_m_upper_limit = photosphere_m_upper_limit
 s% photosphere_m_lower_limit = photosphere_m_lower_limit
 s% photosphere_m_sub_M_center_limit = photosphere_m_sub_M_center_limit
 s% photosphere_r_upper_limit = photosphere_r_upper_limit
 s% photosphere_r_lower_limit = photosphere_r_lower_limit
 s% log_Teff_upper_limit = log_Teff_upper_limit
 s% log_Teff_lower_limit = log_Teff_lower_limit
 s% log_Tsurf_upper_limit = log_Tsurf_upper_limit
 s% log_Tsurf_lower_limit = log_Tsurf_lower_limit
 s% log_Rsurf_upper_limit = log_Rsurf_upper_limit
 s% log_Rsurf_lower_limit = log_Rsurf_lower_limit
 s% log_Psurf_upper_limit = log_Psurf_upper_limit
 s% log_Psurf_lower_limit = log_Psurf_lower_limit
 s% log_Dsurf_upper_limit = log_Dsurf_upper_limit
 s% log_Dsurf_lower_limit = log_Dsurf_lower_limit
 s% log_L_upper_limit = log_L_upper_limit
 s% log_L_lower_limit = log_L_lower_limit
 s% log_g_upper_limit = log_g_upper_limit
 s% log_g_lower_limit = log_g_lower_limit

 s% power_nuc_burn_upper_limit = power_nuc_burn_upper_limit
 s% power_h_burn_upper_limit = power_h_burn_upper_limit
 s% power_he_burn_upper_limit = power_he_burn_upper_limit
 s% power_z_burn_upper_limit = power_z_burn_upper_limit
 s% power_nuc_burn_lower_limit = power_nuc_burn_lower_limit
 s% power_h_burn_lower_limit = power_h_burn_lower_limit
 s% power_he_burn_lower_limit = power_he_burn_lower_limit
 s% power_z_burn_lower_limit = power_z_burn_lower_limit

 ! output of "snapshots" for restarts
 s% photo_interval = photo_interval
 s% photo_digits = photo_digits
 s% photo_directory = photo_directory
 ! output of history and profiles.
 s% do_history_file = do_history_file
 s% history_interval = history_interval

 s% write_header_frequency = write_header_frequency
 s% terminal_interval = terminal_interval
 s% terminal_show_age_units = terminal_show_age_units
 s% terminal_show_timestep_units = terminal_show_timestep_units
 s% terminal_show_log_dt = terminal_show_log_dt
 s% terminal_show_log_age = terminal_show_log_age
 s% extra_terminal_output_file = extra_terminal_output_file
 s% num_trace_history_values = num_trace_history_values
 s% trace_history_value_name = trace_history_value_name

 s% log_directory = log_directory
 s% star_history_name = star_history_name
 s% star_history_header_name = star_history_header_name
 s% star_history_dbl_format = star_history_dbl_format
 s% star_history_int_format = star_history_int_format
 s% star_history_txt_format = star_history_txt_format

 s% profiles_index_name = profiles_index_name
 s% profile_data_prefix = profile_data_prefix
 s% profile_data_suffix = profile_data_suffix
 s% profile_data_header_suffix = profile_data_header_suffix
 s% profile_int_format = profile_int_format
 s% profile_txt_format = profile_txt_format
 s% profile_dbl_format = profile_dbl_format
 s% profile_header_include_sys_details = profile_header_include_sys_details
 s% write_profiles_flag = write_profiles_flag
 s% profile_interval = profile_interval
 s% priority_profile_interval = priority_profile_interval
 s% profile_model = profile_model
 s% max_num_profile_models = max_num_profile_models
 s% max_num_profile_zones = max_num_profile_zones

 s% write_controls_info_with_profile = write_controls_info_with_profile
 s% controls_data_prefix = controls_data_prefix
 s% controls_data_suffix = controls_data_suffix

 s% write_pulse_data_with_profile = write_pulse_data_with_profile
 s% pulse_data_format = pulse_data_format
 s% add_atmosphere_to_pulse_data = add_atmosphere_to_pulse_data
 s% add_center_point_to_pulse_data = add_center_point_to_pulse_data
 s% keep_surface_point_for_pulse_data = keep_surface_point_for_pulse_data
 s% add_double_points_to_pulse_data = add_double_points_to_pulse_data
 s% interpolate_rho_for_pulse_data = interpolate_rho_for_pulse_data
 s% threshold_grad_mu_for_double_point = threshold_grad_mu_for_double_point
 s% max_number_of_double_points = max_number_of_double_points

 s% fgong_header = fgong_header
 s% fgong_ivers = fgong_ivers

 s% max_num_gyre_points = max_num_gyre_points
 s% format_for_OSC_data = format_for_OSC_data
 s% fgong_zero_A_inside_r = fgong_zero_A_inside_r
 s% use_other_export_pulse_data = use_other_export_pulse_data
 s% use_other_get_pulse_data = use_other_get_pulse_data
 s% use_other_edit_pulse_data = use_other_edit_pulse_data

 s% write_model_with_profile = write_model_with_profile
 s% model_data_prefix = model_data_prefix
 s% model_data_suffix = model_data_suffix

 s% mixing_D_limit_for_log = mixing_D_limit_for_log
 s% trace_mass_location = trace_mass_location
 s% min_tau_for_max_abs_v_location = min_tau_for_max_abs_v_location
 s% min_q_for_inner_mach1_location = min_q_for_inner_mach1_location
 s% max_q_for_outer_mach1_location = max_q_for_outer_mach1_location
 
 s% mass_depth_for_L_surf = mass_depth_for_L_surf
 s% conv_core_gap_dq_limit = conv_core_gap_dq_limit

 ! burn zone eps definitions for use in logs and profiles
 s% burn_min1 = burn_min1
 s% burn_min2 = burn_min2

 s% max_conv_vel_div_csound_maxq = max_conv_vel_div_csound_maxq
 s% width_for_limit_conv_vel = width_for_limit_conv_vel
 s% max_q_for_limit_conv_vel = max_q_for_limit_conv_vel
 s% max_mass_in_gm_for_limit_conv_vel = max_mass_in_gm_for_limit_conv_vel
 s% max_r_in_cm_for_limit_conv_vel = max_r_in_cm_for_limit_conv_vel

 ! for reported average values
 s% surface_avg_abundance_dq = surface_avg_abundance_dq
 s% center_avg_value_dq = center_avg_value_dq

 ! mixing parameters
 s% min_convective_gap = min_convective_gap
 s% min_thermohaline_gap = min_thermohaline_gap
 s% min_semiconvection_gap = min_semiconvection_gap
 s% min_thermohaline_dropout = min_thermohaline_dropout
 s% max_dropout_gradL_sub_grada = max_dropout_gradL_sub_grada
 s% remove_embedded_semiconvection = remove_embedded_semiconvection
 s% recalc_mix_info_after_evolve = recalc_mix_info_after_evolve
 s% remove_mixing_glitches = remove_mixing_glitches
 s% okay_to_remove_mixing_singleton = okay_to_remove_mixing_singleton
 s% prune_bad_cz_min_Hp_height = prune_bad_cz_min_Hp_height
 s% prune_bad_cz_min_log_eps_nuc = prune_bad_cz_min_log_eps_nuc
 s% redo_conv_for_dr_lt_mixing_length = redo_conv_for_dr_lt_mixing_length

 s% alpha_semiconvection = alpha_semiconvection
 s% semiconvection_option = semiconvection_option
 s% use_Ledoux_criterion = use_Ledoux_criterion
 s% num_cells_for_smooth_gradL_composition_term = num_cells_for_smooth_gradL_composition_term
 s% threshold_for_smooth_gradL_composition_term = threshold_for_smooth_gradL_composition_term
 s% clip_D_limit = clip_D_limit
 s% fix_eps_grav_transition_to_grid = fix_eps_grav_transition_to_grid

s% okay_to_reduce_gradT_excess = okay_to_reduce_gradT_excess
s% gradT_excess_f1 = gradT_excess_f1
s% gradT_excess_f2 = gradT_excess_f2
s% gradT_excess_age_fraction = gradT_excess_age_fraction
s% gradT_excess_max_change = gradT_excess_max_change
s% gradT_excess_lambda1 = gradT_excess_lambda1
s% gradT_excess_beta1 = gradT_excess_beta1
s% gradT_excess_lambda2 = gradT_excess_lambda2
s% gradT_excess_beta2 = gradT_excess_beta2
s% gradT_excess_dlambda = gradT_excess_dlambda
s% gradT_excess_dbeta = gradT_excess_dbeta
s% gradT_excess_max_center_h1 = gradT_excess_max_center_h1
s% gradT_excess_min_center_he4 = gradT_excess_min_center_he4
s% gradT_excess_max_logT = gradT_excess_max_logT
s% gradT_excess_min_log_tau_full_on = gradT_excess_min_log_tau_full_on
s% gradT_excess_max_log_tau_full_off = gradT_excess_max_log_tau_full_off

 s% D_mix_zero_region_bottom_q = D_mix_zero_region_bottom_q
 s% D_mix_zero_region_top_q = D_mix_zero_region_top_q
 s% dq_D_mix_zero_at_H_He_crossover = dq_D_mix_zero_at_H_He_crossover
 s% dq_D_mix_zero_at_H_C_crossover = dq_D_mix_zero_at_H_C_crossover

 s% use_superad_reduction = use_superad_reduction
 s% superad_reduction_gamma_limit = superad_reduction_gamma_limit
 s% superad_reduction_gamma_limit_scale = superad_reduction_gamma_limit_scale
 s% superad_reduction_gamma_inv_scale = superad_reduction_gamma_inv_scale
 s% superad_reduction_diff_grads_limit = superad_reduction_diff_grads_limit
 s% superad_reduction_limit = superad_reduction_limit
 
 s% max_logT_for_mlt = max_logT_for_mlt
 s% mlt_make_surface_no_mixing = mlt_make_surface_no_mixing
 s% do_normalize_dqs_as_part_of_set_qs = do_normalize_dqs_as_part_of_set_qs

 s% thermohaline_coeff = thermohaline_coeff
 s% thermohaline_option = thermohaline_option
 s% mixing_length_alpha = mixing_length_alpha
 s% remove_small_D_limit = remove_small_D_limit
 s% alt_scale_height_flag = alt_scale_height_flag
 s% Henyey_MLT_y_param = Henyey_MLT_y_param
 s% Henyey_MLT_nu_param = Henyey_MLT_nu_param
 s% make_gradr_sticky_in_solver_iters = make_gradr_sticky_in_solver_iters
 s% min_logT_for_make_gradr_sticky_in_solver_iters = min_logT_for_make_gradr_sticky_in_solver_iters
 s% no_MLT_below_shock = no_MLT_below_shock
 s% MLT_option = MLT_option
 s% mlt_use_rotation_correction = mlt_use_rotation_correction
 s% mlt_Pturb_factor = mlt_Pturb_factor

 s% burn_z_mix_region_logT = burn_z_mix_region_logT
 s% burn_he_mix_region_logT = burn_he_mix_region_logT
 s% burn_h_mix_region_logT = burn_h_mix_region_logT
 s% max_Y_for_burn_z_mix_region = max_Y_for_burn_z_mix_region
 s% max_X_for_burn_he_mix_region = max_X_for_burn_he_mix_region
 
 s% limit_overshoot_Hp_using_size_of_convection_zone = limit_overshoot_Hp_using_size_of_convection_zone

 s%predictive_mix = predictive_mix
 s%predictive_superad_thresh = predictive_superad_thresh
 s%predictive_avoid_reversal = predictive_avoid_reversal
 s%predictive_limit_ingestion = predictive_limit_ingestion
 s%predictive_ingestion_factor = predictive_ingestion_factor
 s%predictive_zone_type = predictive_zone_type
 s%predictive_zone_loc = predictive_zone_loc
 s%predictive_bdy_loc = predictive_bdy_loc
 s%predictive_bdy_q_min = predictive_bdy_q_min
 s%predictive_bdy_q_max = predictive_bdy_q_max

 s%do_conv_premix = do_conv_premix
 s%conv_premix_avoid_increase = conv_premix_avoid_increase
 s%conv_premix_time_factor = conv_premix_time_factor
 s%conv_premix_fix_pgas = conv_premix_fix_pgas
 s%conv_premix_dump_snapshots = conv_premix_dump_snapshots
 s%do_premix_heating = do_premix_heating

 s%overshoot_f = overshoot_f
 s%overshoot_f0 = overshoot_f0
 s%overshoot_D0 = overshoot_D0
 s%overshoot_Delta0 = overshoot_Delta0
 s%overshoot_mass_full_on = overshoot_mass_full_on
 s%overshoot_mass_full_off = overshoot_mass_full_off
 s%overshoot_scheme = overshoot_scheme
 s%overshoot_zone_type = overshoot_zone_type
 s%overshoot_zone_loc = overshoot_zone_loc
 s%overshoot_bdy_loc = overshoot_bdy_loc
 s%overshoot_D_min = overshoot_D_min
 s%overshoot_brunt_B_max = overshoot_brunt_B_max

 s% max_conv_vel_div_csound = max_conv_vel_div_csound
 s% max_v_for_convection = max_v_for_convection
 s% max_q_for_convection_with_hydro_on = max_q_for_convection_with_hydro_on
 s% max_v_div_cs_for_convection = max_v_div_cs_for_convection
 s% max_abs_du_div_cs_for_convection = max_abs_du_div_cs_for_convection

 s% calculate_Brunt_B = calculate_Brunt_B
 s% calculate_Brunt_N2 = calculate_Brunt_N2
 s% brunt_N2_coefficient = brunt_N2_coefficient
 s% num_cells_for_smooth_brunt_B = num_cells_for_smooth_brunt_B
 s% threshold_for_smooth_brunt_B = threshold_for_smooth_brunt_B
 s% min_magnitude_brunt_B = min_magnitude_brunt_B

 s% min_overshoot_q = min_overshoot_q
 s% overshoot_alpha = overshoot_alpha

   s% RSP_max_num_periods = RSP_max_num_periods
   s% RSP_target_steps_per_cycle = RSP_target_steps_per_cycle
   s% RSP_min_max_R_for_periods = RSP_min_max_R_for_periods
   s% RSP_min_deltaR_for_periods = RSP_min_deltaR_for_periods
   s% RSP_default_PERIODLIN = RSP_default_PERIODLIN
   s% RSP_min_PERIOD_div_PERIODLIN = RSP_min_PERIOD_div_PERIODLIN
   s% RSP_GREKM_avg_abs_frac_new = RSP_GREKM_avg_abs_frac_new
   s% RSP_GREKM_avg_abs_limit = RSP_GREKM_avg_abs_limit
   s% RSP_theta = RSP_theta
   s% RSP_thetat = RSP_thetat
   s% RSP_thetau = RSP_thetau
   s% RSP_thetae = RSP_thetae
   s% RSP_thetaq = RSP_thetaq
   s% RSP_wtr = RSP_wtr
   s% RSP_wtc = RSP_wtc
   s% RSP_wtt = RSP_wtt
   s% RSP_gam = RSP_gam
   s% RSP_alfa = RSP_alfa
   s% RSP_alfap = RSP_alfap
   s% RSP_alfam = RSP_alfam
   s% RSP_alfat = RSP_alfat
   s% RSP_alfas = RSP_alfas
   s% RSP_alfac = RSP_alfac
   s% RSP_alfad =  RSP_alfad
   s% RSP_gammar = RSP_gammar
   s% RSP_efl0 = RSP_efl0
   s% RSP_min_tau_for_turbulent_flux = RSP_min_tau_for_turbulent_flux
   s% RSP_cq = RSP_cq
   s% RSP_zsh = RSP_zsh
   s% RSP_Qvisc_quadratic = RSP_Qvisc_quadratic
   s% RSP_Qvisc_linear = RSP_Qvisc_linear
   s% RSP_Qvisc_linear_static = RSP_Qvisc_linear_static
   s% RSP_tol_max_corr = RSP_tol_max_corr
   s% RSP_tol_max_resid = RSP_tol_max_resid
   s% RSP_max_iters_per_try = RSP_max_iters_per_try
   s% RSP_max_retries_per_step = RSP_max_retries_per_step
   s% RSP_nz_div_IBOTOM = RSP_nz_div_IBOTOM
   s% RSP_kick_vsurf_km_per_sec = RSP_kick_vsurf_km_per_sec
   s% RSP_fraction_1st_overtone = RSP_fraction_1st_overtone
   s% RSP_fraction_2nd_overtone = RSP_fraction_2nd_overtone
   s% RSP_Avel = RSP_Avel
   s% RSP_Arnd = RSP_Arnd
   s% RSP_mode_for_setting_PERIODLIN = RSP_mode_for_setting_PERIODLIN
   s% RSP_initial_dt_factor = RSP_initial_dt_factor
   s% RSP_v_div_cs_threshold_for_dt_limit = RSP_v_div_cs_threshold_for_dt_limit
   s% RSP_max_dt_times_min_dr_div_cs = RSP_max_dt_times_min_dr_div_cs
   s% RSP_max_dt_times_min_rad_diff_time = RSP_max_dt_times_min_rad_diff_time
   s% RSP_max_dt = RSP_max_dt
   s% RSP_testing = RSP_testing
   s% RSP_report_limit_dt = RSP_report_limit_dt
   s% RSP_use_Prad_for_Psurf = RSP_use_Prad_for_Psurf
   s% RSP_report_undercorrections = RSP_report_undercorrections
   s% RSP_use_atm_grey_with_kap_for_Psurf = RSP_use_atm_grey_with_kap_for_Psurf
   s% use_other_RSP_linear_analysis = use_other_RSP_linear_analysis
   s% use_other_RSP_build_model = use_other_RSP_build_model
   s% RSP_kap_density_factor = RSP_kap_density_factor
   s% RSP_fixed_Psurf = RSP_fixed_Psurf
   s% RSP_hydro_only = RSP_hydro_only
   s% RSP_tau_surf_for_atm_grey_with_kap = RSP_tau_surf_for_atm_grey_with_kap
   s% RSP_Psurf = RSP_Psurf
   s% set_RSP_Psurf_to_multiple_of_initial_P1 = set_RSP_Psurf_to_multiple_of_initial_P1
   s% RSP_surface_tau = RSP_surface_tau
   s% RSP_write_map = RSP_write_map
   s% RSP_trace_RSP_build_model = RSP_trace_RSP_build_model
   s% RSP_map_filename = RSP_map_filename
   s% RSP_map_columns_filename = RSP_map_columns_filename
   s% RSP_map_history_filename = RSP_map_history_filename
   s% RSP_map_first_period = RSP_map_first_period
   s% RSP_map_last_period = RSP_map_last_period
   s% RSP_map_zone_interval = RSP_map_zone_interval
   s% RSP_nmodes = RSP_nmodes
   s% RSP_work_period = RSP_work_period
   s% RSP_work_filename = RSP_work_filename
   s% RSP_nz_outer = RSP_nz_outer
   s% RSP_max_outer_dm_tries = RSP_max_outer_dm_tries
   s% RSP_max_inner_scale_tries = RSP_max_inner_scale_tries
   s% RSP_relax_max_tries = RSP_relax_max_tries
   s% RSP_T_anchor_tolerance = RSP_T_anchor_tolerance
   s% RSP_T_inner_tolerance = RSP_T_inner_tolerance
   s% RSP_relax_dm_tolerance = RSP_relax_dm_tolerance
   s% RSP_dq_1_factor = RSP_dq_1_factor
   s% use_RSP_new_start_scheme = use_RSP_new_start_scheme
   s% RSP_do_check_omega = RSP_do_check_omega
   s% RSP_report_check_omega_changes = RSP_report_check_omega_changes
   s% RSP_nz = RSP_nz
   s% RSP_T_anchor = RSP_T_anchor
   s% RSP_T_inner = RSP_T_inner
   s% RSP_relax_initial_model = RSP_relax_initial_model
   s% RSP_relax_alfap_before_alfat = RSP_relax_alfap_before_alfat
   s% RSP_relax_adjust_inner_mass_distribution = RSP_relax_adjust_inner_mass_distribution
   s% RSP_Teff = RSP_Teff
   s% RSP_mass = RSP_mass
   s% RSP_L = RSP_L
   s% RSP_X = RSP_X
   s% RSP_Z = RSP_Z

 s% RTI_smooth_mass = RTI_smooth_mass
 s% RTI_smooth_iterations = RTI_smooth_iterations
 s% RTI_smooth_fraction = RTI_smooth_fraction

 s% alpha_RTI_diffusion_factor = alpha_RTI_diffusion_factor
 s% dudt_RTI_diffusion_factor = dudt_RTI_diffusion_factor
 s% dedt_RTI_diffusion_factor = dedt_RTI_diffusion_factor
 s% dlnddt_RTI_diffusion_factor = dlnddt_RTI_diffusion_factor
 s% composition_RTI_diffusion_factor = composition_RTI_diffusion_factor
 s% max_M_RTI_factors_full_on = max_M_RTI_factors_full_on
 s% min_M_RTI_factors_full_off = min_M_RTI_factors_full_off

 s% alpha_RTI_src_min_v_div_cs = alpha_RTI_src_min_v_div_cs
 s% alpha_RTI_src_max_q = alpha_RTI_src_max_q
 s% alpha_RTI_src_min_q = alpha_RTI_src_min_q

 s% T_mix_limit = T_mix_limit
 s% mlt_gradT_fraction = mlt_gradT_fraction

 ! atmosphere -- surface boundary conditions
 s% atm_option = atm_option
 s% atm_off_table_option = atm_off_table_option
 s% Pextra_factor = Pextra_factor
 s% atm_fixed_Teff = atm_fixed_Teff
 s% atm_fixed_Psurf = atm_fixed_Psurf
 s% atm_fixed_Tsurf = atm_fixed_Tsurf

 s% atm_T_tau_relation = atm_T_tau_relation
 s% atm_T_tau_opacity = atm_T_tau_opacity
 s% atm_T_tau_errtol = atm_T_tau_errtol
 s% atm_T_tau_max_iters = atm_T_tau_max_iters
 s% atm_T_tau_max_steps = atm_T_tau_max_steps

 s% atm_table = atm_table

 s% atm_irradiated_opacity = atm_irradiated_opacity
 s% atm_irradiated_errtol = atm_irradiated_errtol
 s% atm_irradiated_T_eq = atm_irradiated_T_eq
 s% atm_irradiated_kap_v = atm_irradiated_kap_v
 s% atm_irradiated_kap_v_div_kap_th = atm_irradiated_kap_v_div_kap_th
 s% atm_irradiated_P_surf = atm_irradiated_P_surf
 s% atm_irradiated_max_iters = atm_irradiated_max_iters

 s% use_compression_outer_BC = use_compression_outer_BC
 s% use_momentum_outer_BC = use_momentum_outer_BC
 s% Tsurf_factor = Tsurf_factor
 s% use_zero_Pgas_outer_BC = use_zero_Pgas_outer_BC
 s% fixed_vsurf = fixed_vsurf
 s% use_fixed_vsurf_outer_BC = use_fixed_vsurf_outer_BC
 s% fixed_Psurf = fixed_Psurf
 s% use_fixed_Psurf_outer_BC = use_fixed_Psurf_outer_BC

 s% atm_build_tau_outer = atm_build_tau_outer
 s% atm_build_dlogtau = atm_build_dlogtau
 s% atm_build_errtol = atm_build_errtol

 s% use_T_tau_gradr_factor = use_T_tau_gradr_factor

 ! extra heat near surface to model irradiation
 s% irradiation_flux = irradiation_flux
 s% column_depth_for_irradiation = column_depth_for_irradiation

 ! extra heat
 s% inject_uniform_extra_heat = inject_uniform_extra_heat
 s% min_q_for_uniform_extra_heat = min_q_for_uniform_extra_heat
 s% max_q_for_uniform_extra_heat = max_q_for_uniform_extra_heat
 s% inject_extra_ergs_sec = inject_extra_ergs_sec
 s% base_of_inject_extra_ergs_sec = base_of_inject_extra_ergs_sec
 s% total_mass_for_inject_extra_ergs_sec = total_mass_for_inject_extra_ergs_sec
 s% start_time_for_inject_extra_ergs_sec = start_time_for_inject_extra_ergs_sec
 s% duration_for_inject_extra_ergs_sec = duration_for_inject_extra_ergs_sec
 s% inject_until_reach_model_with_total_energy = inject_until_reach_model_with_total_energy

 ! mass gain or loss
 s% mass_change = mass_change
 s% mass_change_full_off_dt = mass_change_full_off_dt
 s% mass_change_full_on_dt = mass_change_full_on_dt
 s% trace_dt_control_mass_change = trace_dt_control_mass_change
 s% no_wind_if_no_rotation = no_wind_if_no_rotation

 s% min_wind = min_wind
 s% max_wind = max_wind
 s% use_accreted_material_j = use_accreted_material_j
 s% accreted_material_j = accreted_material_j
 s% D_omega_mixing_rate = D_omega_mixing_rate
 s% D_omega_mixing_across_convection_boundary = D_omega_mixing_across_convection_boundary
 s% max_q_for_D_omega_zero_in_convection_region = max_q_for_D_omega_zero_in_convection_region
 s% nu_omega_mixing_rate = nu_omega_mixing_rate
 s% nu_omega_mixing_across_convection_boundary = nu_omega_mixing_across_convection_boundary
 s% max_q_for_nu_omega_zero_in_convection_region = max_q_for_nu_omega_zero_in_convection_region

 s% mdot_omega_power = mdot_omega_power
 s% max_rotational_mdot_boost = max_rotational_mdot_boost
 s% max_mdot_jump_for_rotation = max_mdot_jump_for_rotation
 s% lim_trace_rotational_mdot_boost = lim_trace_rotational_mdot_boost
 s% rotational_mdot_boost_fac = rotational_mdot_boost_fac
 s% rotational_mdot_kh_fac = rotational_mdot_kh_fac
 s% surf_avg_tau = surf_avg_tau
 s% surf_avg_tau_min = surf_avg_tau_min

 s% super_eddington_scaling_factor = super_eddington_scaling_factor
 s% super_eddington_wind_Ledd_factor = super_eddington_wind_Ledd_factor
 s% wind_boost_full_off_L_div_Ledd = wind_boost_full_off_L_div_Ledd
 s% wind_boost_full_on_L_div_Ledd = wind_boost_full_on_L_div_Ledd
 s% super_eddington_wind_max_boost = super_eddington_wind_max_boost
 s% trace_super_eddington_wind_boost = trace_super_eddington_wind_boost
 
 s% max_tries_for_implicit_wind = max_tries_for_implicit_wind
 s% iwind_tolerance = iwind_tolerance
 s% iwind_lambda = iwind_lambda

 s% cool_wind_full_on_T = cool_wind_full_on_T
 s% hot_wind_full_on_T = hot_wind_full_on_T

 s% rlo_scaling_factor = rlo_scaling_factor
 s% rlo_wind_min_L = rlo_wind_min_L
 s% rlo_wind_max_Teff = rlo_wind_max_Teff
 s% rlo_wind_roche_lobe_radius = rlo_wind_roche_lobe_radius
 s% roche_lobe_xfer_full_on = roche_lobe_xfer_full_on
 s% roche_lobe_xfer_full_off = roche_lobe_xfer_full_off
 s% rlo_wind_base_mdot = rlo_wind_base_mdot
 s% rlo_wind_scale_height = rlo_wind_scale_height

 s% hot_wind_scheme = hot_wind_scheme
 s% cool_wind_RGB_scheme = cool_wind_RGB_scheme
 s% cool_wind_AGB_scheme = cool_wind_AGB_scheme
 s% RGB_to_AGB_wind_switch = RGB_to_AGB_wind_switch
 s% Reimers_scaling_factor = Reimers_scaling_factor
 s% Blocker_scaling_factor = Blocker_scaling_factor
 s% de_Jager_scaling_factor = de_Jager_scaling_factor
 s% van_Loon_scaling_factor = van_Loon_scaling_factor
 s% Nieuwenhuijzen_scaling_factor = Nieuwenhuijzen_scaling_factor
 s% Vink_scaling_factor = Vink_scaling_factor
 s% Dutch_scaling_factor = Dutch_scaling_factor
 s% Dutch_wind_lowT_scheme = Dutch_wind_lowT_scheme

 s% wind_H_envelope_limit = wind_H_envelope_limit
 s% wind_H_He_envelope_limit = wind_H_He_envelope_limit
 s% wind_He_layer_limit = wind_He_layer_limit

 s% max_logT_for_k_below_const_q = max_logT_for_k_below_const_q
 s% max_q_for_k_below_const_q = max_q_for_k_below_const_q
 s% min_q_for_k_below_const_q = min_q_for_k_below_const_q
 s% max_logT_for_k_const_mass = max_logT_for_k_const_mass
 s% min_q_for_k_const_mass = min_q_for_k_const_mass
 s% max_q_for_k_const_mass = max_q_for_k_const_mass

 ! composition of added mass
 s% accrete_same_as_surface = accrete_same_as_surface

 s% accrete_given_mass_fractions = accrete_given_mass_fractions
 s% num_accretion_species = num_accretion_species
 s% accretion_species_id = accretion_species_id
 s% accretion_species_xa = accretion_species_xa

 s% accretion_h1 = accretion_h1
 s% accretion_h2 = accretion_h2
 s% accretion_he3 = accretion_he3
 s% accretion_he4 = accretion_he4
 s% accretion_zfracs = accretion_zfracs
 s% accretion_dump_missing_metals_into_heaviest = accretion_dump_missing_metals_into_heaviest

 ! special list of z fractions
 s% z_fraction_li = z_fraction_li
 s% z_fraction_be = z_fraction_be
 s% z_fraction_b = z_fraction_b
 s% z_fraction_c = z_fraction_c
 s% z_fraction_n = z_fraction_n
 s% z_fraction_o = z_fraction_o
 s% z_fraction_f = z_fraction_f
 s% z_fraction_ne = z_fraction_ne
 s% z_fraction_na = z_fraction_na
 s% z_fraction_mg = z_fraction_mg
 s% z_fraction_al = z_fraction_al
 s% z_fraction_si = z_fraction_si
 s% z_fraction_p = z_fraction_p
 s% z_fraction_s = z_fraction_s
 s% z_fraction_cl = z_fraction_cl
 s% z_fraction_ar = z_fraction_ar
 s% z_fraction_k = z_fraction_k
 s% z_fraction_ca = z_fraction_ca
 s% z_fraction_sc = z_fraction_sc
 s% z_fraction_ti = z_fraction_ti
 s% z_fraction_v = z_fraction_v
 s% z_fraction_cr = z_fraction_cr
 s% z_fraction_mn = z_fraction_mn
 s% z_fraction_fe = z_fraction_fe
 s% z_fraction_co = z_fraction_co
 s% z_fraction_ni = z_fraction_ni
 s% z_fraction_cu = z_fraction_cu
 s% z_fraction_zn = z_fraction_zn

 s% lgT_lo_for_set_new_abundances = lgT_lo_for_set_new_abundances
 s% lgT_hi_for_set_new_abundances = lgT_hi_for_set_new_abundances

 ! automatic stops for mass loss/gain
 s% max_star_mass_for_gain = max_star_mass_for_gain
 s% min_star_mass_for_loss = min_star_mass_for_loss
 s% max_T_center_for_any_mass_loss = max_T_center_for_any_mass_loss
 s% max_T_center_for_full_mass_loss = max_T_center_for_full_mass_loss

 ! extra power source
 s% extra_power_source = extra_power_source

 ! relaxation parameters
 s% relax_dlnZ = relax_dlnZ
 s% relax_dY = relax_dY

 ! mesh adjustment
 s% show_mesh_changes = show_mesh_changes
 s% okay_to_remesh = okay_to_remesh
 s% restore_mesh_on_retry = restore_mesh_on_retry
 s% num_steps_to_hold_mesh_after_retry = num_steps_to_hold_mesh_after_retry
 s% trace_mesh_adjust_error_in_conservation = trace_mesh_adjust_error_in_conservation
 s% max_rel_delta_IE_for_mesh_total_energy_balance = max_rel_delta_IE_for_mesh_total_energy_balance
 s% max_allowed_nz = max_allowed_nz
 s% mesh_max_allowed_ratio = mesh_max_allowed_ratio
 s% remesh_max_allowed_logT = remesh_max_allowed_logT
 s% max_delta_x_for_merge = max_delta_x_for_merge

 s% mesh_ok_to_merge = mesh_ok_to_merge
 s% mesh_max_k_old_for_split = mesh_max_k_old_for_split
 s% mesh_min_k_old_for_split = mesh_min_k_old_for_split
 s% mesh_adjust_get_T_from_E = mesh_adjust_get_T_from_E

 s% max_dq = max_dq
 s% min_dq = min_dq
 s% min_dq_for_split = min_dq_for_split
 s% min_dq_for_xa = min_dq_for_xa
 s% min_dq_for_xa_convective = min_dq_for_xa_convective
 s% min_dq_for_logT = min_dq_for_logT

 s% mesh_min_dlnR = mesh_min_dlnR
 s% merge_if_dlnR_too_small = merge_if_dlnR_too_small

 s% mesh_min_dr_div_dRstar = mesh_min_dr_div_dRstar
 s% merge_if_dr_div_dRstar_too_small = merge_if_dr_div_dRstar_too_small

 s% mesh_min_dr_div_cs = mesh_min_dr_div_cs
 s% merge_if_dr_div_cs_too_small = merge_if_dr_div_cs_too_small

 s% max_center_cell_dq = max_center_cell_dq
 s% max_surface_cell_dq = max_surface_cell_dq
 s% max_num_subcells = max_num_subcells
 s% max_num_merge_cells = max_num_merge_cells

 s% mesh_delta_coeff = mesh_delta_coeff
 s% mesh_delta_coeff_for_highT = mesh_delta_coeff_for_highT
 s% logT_max_for_standard_mesh_delta_coeff = logT_max_for_standard_mesh_delta_coeff
 s% logT_min_for_highT_mesh_delta_coeff = logT_min_for_highT_mesh_delta_coeff
 s% mesh_Pgas_div_P_exponent = mesh_Pgas_div_P_exponent

 s% remesh_dt_limit = remesh_dt_limit

 s% E_function_weight = E_function_weight
 s% E_function_param = E_function_param
 s% P_function_weight = P_function_weight

 s% mesh_logX_species = mesh_logX_species
 s% mesh_logX_min_for_extra = mesh_logX_min_for_extra
 s% mesh_dlogX_dlogP_extra = mesh_dlogX_dlogP_extra
 s% mesh_dlogX_dlogP_full_on = mesh_dlogX_dlogP_full_on
 s% mesh_dlogX_dlogP_full_off = mesh_dlogX_dlogP_full_off

 s% mesh_dlog_eps_min_for_extra = mesh_dlog_eps_min_for_extra
 s% mesh_dlog_eps_dlogP_full_on = mesh_dlog_eps_dlogP_full_on
 s% mesh_dlog_eps_dlogP_full_off = mesh_dlog_eps_dlogP_full_off

 s% mesh_dlog_pp_dlogP_extra = mesh_dlog_pp_dlogP_extra
 s% mesh_dlog_cno_dlogP_extra = mesh_dlog_cno_dlogP_extra
 s% mesh_dlog_3alf_dlogP_extra = mesh_dlog_3alf_dlogP_extra

 s% mesh_dlog_burn_c_dlogP_extra = mesh_dlog_burn_c_dlogP_extra
 s% mesh_dlog_burn_n_dlogP_extra = mesh_dlog_burn_n_dlogP_extra
 s% mesh_dlog_burn_o_dlogP_extra = mesh_dlog_burn_o_dlogP_extra
 s% mesh_dlog_burn_ne_dlogP_extra = mesh_dlog_burn_ne_dlogP_extra
 s% mesh_dlog_burn_na_dlogP_extra = mesh_dlog_burn_na_dlogP_extra
 s% mesh_dlog_burn_mg_dlogP_extra = mesh_dlog_burn_mg_dlogP_extra
 s% mesh_dlog_burn_si_dlogP_extra = mesh_dlog_burn_si_dlogP_extra
 s% mesh_dlog_burn_s_dlogP_extra = mesh_dlog_burn_s_dlogP_extra
 s% mesh_dlog_burn_ar_dlogP_extra = mesh_dlog_burn_ar_dlogP_extra
 s% mesh_dlog_burn_ca_dlogP_extra = mesh_dlog_burn_ca_dlogP_extra
 s% mesh_dlog_burn_ti_dlogP_extra = mesh_dlog_burn_ti_dlogP_extra
 s% mesh_dlog_burn_cr_dlogP_extra = mesh_dlog_burn_cr_dlogP_extra
 s% mesh_dlog_burn_fe_dlogP_extra = mesh_dlog_burn_fe_dlogP_extra

 s% mesh_dlog_cc_dlogP_extra = mesh_dlog_cc_dlogP_extra
 s% mesh_dlog_co_dlogP_extra = mesh_dlog_co_dlogP_extra
 s% mesh_dlog_oo_dlogP_extra = mesh_dlog_oo_dlogP_extra

 s% mesh_dlog_pnhe4_dlogP_extra = mesh_dlog_pnhe4_dlogP_extra
 s% mesh_dlog_photo_dlogP_extra = mesh_dlog_photo_dlogP_extra
 s% mesh_dlog_other_dlogP_extra = mesh_dlog_other_dlogP_extra
 
 s% mesh_delta_coeff_factor_smooth_iters = mesh_delta_coeff_factor_smooth_iters

 s% T_function1_weight = T_function1_weight
 s% T_function2_weight = T_function2_weight
 s% T_function2_param = T_function2_param

 s% R_function_weight = R_function_weight
 s% R_function_param = R_function_param

 s% R_function2_weight = R_function2_weight
 s% R_function2_param1 = R_function2_param1
 s% R_function2_param2 = R_function2_param2

 s% R_function3_weight = R_function3_weight

 s% M_function_weight = M_function_weight
 s% M_function_param = M_function_param

 s% gradT_function_weight = gradT_function_weight
 s% log_tau_function_weight = log_tau_function_weight
 s% log_kap_function_weight = log_kap_function_weight
 s% omega_function_weight = omega_function_weight

 s% gam_function_weight = gam_function_weight
 s% gam_function_param1 = gam_function_param1
 s% gam_function_param2 = gam_function_param2

 s% xa_function_species = xa_function_species
 s% xa_function_weight = xa_function_weight
 s% xa_function_param = xa_function_param
 s% xa_mesh_delta_coeff = xa_mesh_delta_coeff
 
 s% use_split_merge_amr = use_split_merge_amr
 s% split_merge_amr_nz_baseline = split_merge_amr_nz_baseline
 s% split_merge_amr_nz_r_core = split_merge_amr_nz_r_core
 s% split_merge_amr_nz_r_core_fraction = split_merge_amr_nz_r_core_fraction
 s% split_merge_amr_mesh_delta_coeff = split_merge_amr_mesh_delta_coeff
 s% split_merge_amr_log_zoning = split_merge_amr_log_zoning
 s% split_merge_amr_hybrid_zoning = split_merge_amr_hybrid_zoning
 s% split_merge_amr_flipped_hybrid_zoning = split_merge_amr_flipped_hybrid_zoning
 s% split_merge_amr_logtau_zoning = split_merge_amr_logtau_zoning
 s% split_merge_amr_okay_to_split_nz = split_merge_amr_okay_to_split_nz
 s% split_merge_amr_okay_to_split_1 = split_merge_amr_okay_to_split_1
 s% merge_amr_inhibit_at_jumps = merge_amr_inhibit_at_jumps
 s% split_merge_amr_MaxLong = split_merge_amr_MaxLong
 s% split_merge_amr_MaxShort = split_merge_amr_MaxShort
 s% merge_amr_max_abs_du_div_cs = merge_amr_max_abs_du_div_cs
 s% merge_amr_ignore_surface_cells = merge_amr_ignore_surface_cells
 s% merge_amr_du_div_cs_limit_only_for_compression = merge_amr_du_div_cs_limit_only_for_compression
 s% split_merge_amr_avoid_repeated_remesh = split_merge_amr_avoid_repeated_remesh
 s% merge_amr_k_for_ignore_surface_cells = merge_amr_k_for_ignore_surface_cells
 s% split_merge_amr_dq_min = split_merge_amr_dq_min
 s% split_merge_amr_dq_max = split_merge_amr_dq_max
 s% split_merge_amr_r_core_cm = split_merge_amr_r_core_cm
 s% split_merge_amr_max_iters = split_merge_amr_max_iters
 s% trace_split_merge_amr = trace_split_merge_amr
 s% equal_split_density_amr = equal_split_density_amr

 ! nuclear reaction parameters
 s% screening_mode = screening_mode
 s% default_net_name = default_net_name

 s% net_logTcut_lo = net_logTcut_lo
 s% net_logTcut_lim = net_logTcut_lim

 s% eps_nuc_factor = eps_nuc_factor
 s% op_split_burn_eps_nuc_infall_limit = op_split_burn_eps_nuc_infall_limit
 s% eps_WD_sedimentation_factor = eps_WD_sedimentation_factor
 s% max_abs_eps_nuc = max_abs_eps_nuc
 s% dxdt_nuc_factor = dxdt_nuc_factor
 s% max_abar_for_burning = max_abar_for_burning
 s% fe56ec_fake_factor = fe56ec_fake_factor
 s% min_T_for_fe56ec_fake_factor = min_T_for_fe56ec_fake_factor
 s% weak_rate_factor = weak_rate_factor

 ! mixing
 s% mix_factor = mix_factor

 s% sig_term_limit = sig_term_limit

 s% sig_min_factor_for_high_Tcenter = sig_min_factor_for_high_Tcenter
 s% Tcenter_min_for_sig_min_factor_full_on = Tcenter_min_for_sig_min_factor_full_on
 s% Tcenter_max_for_sig_min_factor_full_off = Tcenter_max_for_sig_min_factor_full_off
 s% max_delta_m_to_bdy_for_sig_min_factor = max_delta_m_to_bdy_for_sig_min_factor
 s% delta_m_lower_for_sig_min_factor = delta_m_lower_for_sig_min_factor
 s% delta_m_upper_for_sig_min_factor = delta_m_upper_for_sig_min_factor

 s% am_sig_term_limit = am_sig_term_limit
 s% am_D_mix_factor = am_D_mix_factor
 s% am_gradmu_factor = am_gradmu_factor
 s% am_nu_factor = am_nu_factor

 s% D_visc_factor = D_visc_factor
 s% D_DSI_factor = D_DSI_factor
 s% D_SH_factor = D_SH_factor
 s% D_SSI_factor = D_SSI_factor
 s% D_ES_factor = D_ES_factor
 s% D_GSF_factor = D_GSF_factor
 s% D_ST_factor = D_ST_factor

 s% am_nu_non_rotation_factor = am_nu_non_rotation_factor
 s% skip_rotation_in_convection_zones = skip_rotation_in_convection_zones
 s% am_nu_DSI_factor = am_nu_DSI_factor
 s% am_nu_SH_factor = am_nu_SH_factor
 s% am_nu_SSI_factor = am_nu_SSI_factor
 s% am_nu_ES_factor = am_nu_ES_factor
 s% am_nu_GSF_factor = am_nu_GSF_factor
 s% am_nu_ST_factor = am_nu_ST_factor
 s% am_nu_visc_factor = am_nu_visc_factor

 s% am_nu_omega_rot_factor = am_nu_omega_rot_factor
 s% am_nu_omega_non_rot_factor = am_nu_omega_non_rot_factor
 s% am_nu_j_rot_factor = am_nu_j_rot_factor
 s% am_nu_j_non_rot_factor = am_nu_j_non_rot_factor

 s% smooth_nu_ST = smooth_nu_ST
 s% smooth_D_ST = smooth_D_ST
 s% smooth_D_DSI = smooth_D_DSI
 s% smooth_D_SSI = smooth_D_SSI
 s% smooth_D_SH = smooth_D_SH
 s% smooth_D_GSF = smooth_D_GSF
 s% smooth_D_ES = smooth_D_ES
 s% smooth_D_omega = smooth_D_omega
 s% smooth_am_nu_rot = smooth_am_nu_rot

 s% simple_i_rot_flag = simple_i_rot_flag
 s% do_adjust_J_lost = do_adjust_J_lost
 s% premix_omega = premix_omega
 s% angular_momentum_error_warn = angular_momentum_error_warn
 s% angular_momentum_error_retry = angular_momentum_error_retry
 s% recalc_mixing_info_each_substep = recalc_mixing_info_each_substep
 s% adjust_J_fraction = adjust_J_fraction
 s% min_q_for_adjust_J_lost = min_q_for_adjust_J_lost
 s% min_J_div_delta_J = min_J_div_delta_J
 s% max_mdot_redo_cnt = max_mdot_redo_cnt
 s% mdot_revise_factor = mdot_revise_factor
 s% implicit_mdot_boost = implicit_mdot_boost
 s% min_years_dt_for_redo_mdot = min_years_dt_for_redo_mdot
 s% surf_omega_div_omega_crit_limit = surf_omega_div_omega_crit_limit
 s% surf_omega_div_omega_crit_tol = surf_omega_div_omega_crit_tol
 s% fitted_fp_ft_i_rot = fitted_fp_ft_i_rot 
 s% w_div_wcrit_max = w_div_wcrit_max
 s% w_div_wcrit_max2 = w_div_wcrit_max2
 s% fp_min = fp_min
 s% ft_min = ft_min
 s% fp_error_limit = fp_error_limit
 s% ft_error_limit = ft_error_limit

 s% D_mix_rotation_max_logT_full_on = D_mix_rotation_max_logT_full_on
 s% D_mix_rotation_min_logT_full_off = D_mix_rotation_min_logT_full_off

 s% set_uniform_am_nu_non_rot = set_uniform_am_nu_non_rot
 s% uniform_am_nu_non_rot = uniform_am_nu_non_rot

 s% set_min_am_nu_non_rot = set_min_am_nu_non_rot
 s% min_am_nu_non_rot = min_am_nu_non_rot
 s% min_center_Ye_for_min_am_nu_non_rot = min_center_Ye_for_min_am_nu_non_rot

 s% set_min_D_mix = set_min_D_mix
 s% mass_lower_limit_for_min_D_mix = mass_lower_limit_for_min_D_mix
 s% mass_upper_limit_for_min_D_mix = mass_upper_limit_for_min_D_mix
 s% min_D_mix = min_D_mix
 s% set_min_D_mix_below_Tmax = set_min_D_mix_below_Tmax
 s% min_D_mix_below_Tmax = min_D_mix_below_Tmax
 s% set_min_D_mix_in_H_He = set_min_D_mix_in_H_He
 s% min_D_mix_in_H_He = min_D_mix_in_H_He
 s% min_center_Ye_for_min_D_mix = min_center_Ye_for_min_D_mix
 s% reaction_neuQs_factor = reaction_neuQs_factor
 s% nonlocal_NiCo_kap_gamma = nonlocal_NiCo_kap_gamma
 s% nonlocal_NiCo_decay_heat = nonlocal_NiCo_decay_heat
 s% dtau_gamma_NiCo_decay_heat = dtau_gamma_NiCo_decay_heat
 s% max_logT_for_net = max_logT_for_net
 s% smooth_outer_xa_big = smooth_outer_xa_big
 s% smooth_outer_xa_small = smooth_outer_xa_small

 ! element diffusion parameters
 s% diffusion_use_iben_macdonald = diffusion_use_iben_macdonald
 s% diffusion_use_paquette = diffusion_use_paquette
 s% diffusion_use_cgs_solver = diffusion_use_cgs_solver
 s% diffusion_use_full_net = diffusion_use_full_net
 s% do_WD_sedimentation_heating = do_WD_sedimentation_heating
 s% min_xa_for_WD_sedimentation_heating = min_xa_for_WD_sedimentation_heating
 s% do_diffusion_heating = do_diffusion_heating
 s% do_element_diffusion = do_element_diffusion
 s% cgs_thermal_diffusion_eta_full_on = cgs_thermal_diffusion_eta_full_on
 s% cgs_thermal_diffusion_eta_full_off = cgs_thermal_diffusion_eta_full_off
 s% diffusion_min_dq_at_surface = diffusion_min_dq_at_surface
 s% diffusion_min_T_at_surface = diffusion_min_T_at_surface
 s% diffusion_min_dq_ratio_at_surface = diffusion_min_dq_ratio_at_surface
 s% diffusion_dt_limit = diffusion_dt_limit

 s% diffusion_min_X_hard_limit = diffusion_min_X_hard_limit
 s% diffusion_X_total_atol = diffusion_X_total_atol
 s% diffusion_X_total_rtol = diffusion_X_total_rtol
 s% diffusion_upwind_abs_v_limit = diffusion_upwind_abs_v_limit
 s% diffusion_dt_div_timescale = diffusion_dt_div_timescale
 s% diffusion_min_num_substeps = diffusion_min_num_substeps
 s% diffusion_max_iters_per_substep = diffusion_max_iters_per_substep
 s% diffusion_max_retries_per_substep = diffusion_max_retries_per_substep
 s% diffusion_v_max = diffusion_v_max
 s% diffusion_gamma_full_off = diffusion_gamma_full_off
 s% diffusion_gamma_full_on = diffusion_gamma_full_on
 s% diffusion_T_full_off = diffusion_T_full_off
 s% D_mix_ignore_diffusion = D_mix_ignore_diffusion
 s% diffusion_T_full_on = diffusion_T_full_on
 s% diffusion_calculates_ionization = diffusion_calculates_ionization
 s% diffusion_nsmooth_typical_charge = diffusion_nsmooth_typical_charge
 s% diffusion_tol_correction_max = diffusion_tol_correction_max
 s% diffusion_tol_correction_norm = diffusion_tol_correction_norm

 s% diffusion_AD_dm_full_on = diffusion_AD_dm_full_on
 s% diffusion_AD_dm_full_off = diffusion_AD_dm_full_off
 s% diffusion_AD_boost_factor = diffusion_AD_boost_factor

 s% diffusion_SIG_factor = diffusion_SIG_factor
 s% diffusion_GT_factor = diffusion_GT_factor

 s% diffusion_Vlimit_dm_full_on = diffusion_Vlimit_dm_full_on
 s% diffusion_Vlimit_dm_full_off = diffusion_Vlimit_dm_full_off
 s% diffusion_Vlimit = diffusion_Vlimit

 s% diffusion_max_T_for_radaccel = diffusion_max_T_for_radaccel
 s% diffusion_min_T_for_radaccel = diffusion_min_T_for_radaccel
 s% diffusion_max_Z_for_radaccel = diffusion_max_Z_for_radaccel
 s% diffusion_min_Z_for_radaccel = diffusion_min_Z_for_radaccel
 s% diffusion_screening_for_radaccel = diffusion_screening_for_radaccel
 s% op_mono_data_path = op_mono_data_path
 s% op_mono_data_cache_filename = op_mono_data_cache_filename

 s% show_diffusion_info = show_diffusion_info
 s% show_diffusion_substep_info = show_diffusion_substep_info
 s% show_diffusion_timing = show_diffusion_timing

 s% diffusion_num_classes = diffusion_num_classes
 s% diffusion_class_representative = diffusion_class_representative
 s% diffusion_class_A_max = diffusion_class_A_max
 s% diffusion_class_typical_charge = diffusion_class_typical_charge
 s% diffusion_class_factor = diffusion_class_factor

 s% diffusion_use_isolve = diffusion_use_isolve
 s% diffusion_rtol_for_isolve = diffusion_rtol_for_isolve
 s% diffusion_atol_for_isolve = diffusion_atol_for_isolve
 s% diffusion_maxsteps_for_isolve = diffusion_maxsteps_for_isolve
 s% diffusion_isolve_solver = diffusion_isolve_solver

 ! eos controls
 s% use_fixed_XZ_for_eos = use_fixed_XZ_for_eos
 s% report_eos_settings_at_start_of_run = report_eos_settings_at_start_of_run
 s% fixed_X_for_eos = fixed_X_for_eos
 s% fixed_Z_for_eos = fixed_Z_for_eos
 s% use_d_eos_dxa = use_d_eos_dxa

 ! opacity controls
 s% use_simple_es_for_kap = use_simple_es_for_kap
 s% use_starting_composition_for_kap = use_starting_composition_for_kap

 s% min_kap_for_dPrad_dm_eqn = min_kap_for_dPrad_dm_eqn
 s% low_logT_op_mono_full_off = low_logT_op_mono_full_off
 s% low_logT_op_mono_full_on = low_logT_op_mono_full_on
 s% high_logT_op_mono_full_off = high_logT_op_mono_full_off
 s% high_logT_op_mono_full_on = high_logT_op_mono_full_on
 s% op_mono_min_X_to_include = op_mono_min_X_to_include
 s% use_op_mono_alt_get_kap = use_op_mono_alt_get_kap
  
 s% include_L_in_correction_limits = include_L_in_correction_limits
 s% include_v_in_correction_limits = include_v_in_correction_limits
 s% include_u_in_correction_limits = include_u_in_correction_limits
 s% include_w_in_correction_limits = include_w_in_correction_limits

 ! asteroseismology controls

 s% get_delta_nu_from_scaled_solar = get_delta_nu_from_scaled_solar
 s% nu_max_sun = nu_max_sun
 s% delta_nu_sun = delta_nu_sun
 s% Teff_sun = Teff_sun
 s% delta_Pg_mode_freq = delta_Pg_mode_freq


 ! hydro parameters
 s% opacity_factor = opacity_factor
 s% opacity_max = opacity_max
 s% min_logT_for_opacity_factor_off = min_logT_for_opacity_factor_off
 s% min_logT_for_opacity_factor_on = min_logT_for_opacity_factor_on
 s% max_logT_for_opacity_factor_on = max_logT_for_opacity_factor_on
 s% max_logT_for_opacity_factor_off = max_logT_for_opacity_factor_off

 s% dxdt_nuc_factor = dxdt_nuc_factor
 s% non_nuc_neu_factor = non_nuc_neu_factor
 s% use_dedt_form_of_energy_eqn = use_dedt_form_of_energy_eqn
 s% always_use_dedt_form_of_energy_eqn = always_use_dedt_form_of_energy_eqn
 s% use_time_centered_eps_grav = use_time_centered_eps_grav
 s% no_dedt_form_during_relax = no_dedt_form_during_relax
 s% always_use_eps_grav_form_of_energy_eqn = always_use_eps_grav_form_of_energy_eqn
 s% dedt_eqn_r_scale = dedt_eqn_r_scale
 s% use_mass_corrections = use_mass_corrections
 s% use_gravity_rotation_correction = use_gravity_rotation_correction
 s% eps_grav_factor = eps_grav_factor
 s% eps_mdot_factor = eps_mdot_factor
 s% include_composition_in_eps_grav = include_composition_in_eps_grav
 s% max_abs_rel_change_surf_lnS = max_abs_rel_change_surf_lnS
 s% max_num_surf_revisions = max_num_surf_revisions
 s% Gamma_lnS_eps_grav_full_off = Gamma_lnS_eps_grav_full_off
 s% Gamma_lnS_eps_grav_full_on = Gamma_lnS_eps_grav_full_on

 s% use_dPrad_dm_form_of_T_gradient_eqn = use_dPrad_dm_form_of_T_gradient_eqn
 s% use_gradT_actual_vs_gradT_MLT_for_T_gradient_eqn = use_gradT_actual_vs_gradT_MLT_for_T_gradient_eqn
 s% include_P_in_velocity_time_centering = include_P_in_velocity_time_centering
 s% include_L_in_velocity_time_centering = include_L_in_velocity_time_centering
 s% use_P_d_1_div_rho_form_of_work_when_time_centering_velocity = use_P_d_1_div_rho_form_of_work_when_time_centering_velocity
 s% steps_before_use_velocity_time_centering = steps_before_use_velocity_time_centering
 s% P_theta_for_velocity_time_centering = P_theta_for_velocity_time_centering
 s% L_theta_for_velocity_time_centering = L_theta_for_velocity_time_centering

 s% RTI_A = RTI_A
 s% RTI_B = RTI_B
 s% RTI_C = RTI_C
 s% RTI_D = RTI_D
 s% RTI_max_alpha = RTI_max_alpha
 s% RTI_C_X_factor = RTI_C_X_factor
 s% RTI_C_X0_frac = RTI_C_X0_frac
 s% RTI_dm_for_center_eta_nondecreasing = RTI_dm_for_center_eta_nondecreasing
 s% RTI_min_dm_behind_shock_for_full_on = RTI_min_dm_behind_shock_for_full_on
 s% RTI_energy_floor = RTI_energy_floor
 s% RTI_D_mix_floor = RTI_D_mix_floor
 s% RTI_min_m_for_D_mix_floor = RTI_min_m_for_D_mix_floor
 s% RTI_log_max_boost = RTI_log_max_boost 
 s% RTI_m_full_boost = RTI_m_full_boost
 s% RTI_m_no_boost = RTI_m_no_boost

 s% conv_vel_D = conv_vel_D
 s% conv_vel_siglimit = conv_vel_siglimit
 s% conv_vel_v0 = conv_vel_v0
 s% min_q_for_normal_mlt_gradT_full_off = min_q_for_normal_mlt_gradT_full_off
 s% max_q_for_normal_mlt_gradT_full_on = max_q_for_normal_mlt_gradT_full_on
 s% conv_vel_ignore_thermohaline = conv_vel_ignore_thermohaline
 s% conv_vel_ignore_semiconvection = conv_vel_ignore_semiconvection
 s% conv_vel_fully_lagrangian = conv_vel_fully_lagrangian
 s% conv_vel_include_homologous_term = conv_vel_include_homologous_term
 s% conv_vel_use_mlt_vc_start = conv_vel_use_mlt_vc_start

 s% velocity_logT_lower_bound = velocity_logT_lower_bound
 s% max_dt_yrs_for_velocity_logT_lower_bound = max_dt_yrs_for_velocity_logT_lower_bound
 s% velocity_q_upper_bound = velocity_q_upper_bound

 ! solvers

 s% tol_correction_norm = tol_correction_norm
 s% tol_max_correction = tol_max_correction
 s% correction_xa_limit = correction_xa_limit

 s% tol_correction_high_T_limit = tol_correction_high_T_limit
 s% tol_correction_norm_high_T = tol_correction_norm_high_T
 s% tol_max_correction_high_T = tol_max_correction_high_T

 s% tol_correction_extreme_T_limit = tol_correction_extreme_T_limit
 s% tol_correction_norm_extreme_T = tol_correction_norm_extreme_T
 s% tol_max_correction_extreme_T = tol_max_correction_extreme_T
 
 s% tol_bad_max_correction = tol_bad_max_correction
 s% bad_max_correction_series_limit = bad_max_correction_series_limit

 s% tol_residual_norm1 = tol_residual_norm1
 s% tol_max_residual1 = tol_max_residual1
 s% tol_residual_norm2 = tol_residual_norm2
 s% tol_max_residual2 = tol_max_residual2
 s% tol_residual_norm3 = tol_residual_norm3
 s% tol_max_residual3 = tol_max_residual3
 s% warning_limit_for_max_residual = warning_limit_for_max_residual
 s% trace_solver_damping = trace_solver_damping
 
 s% relax_use_gold_tolerances = relax_use_gold_tolerances
 s% relax_tol_correction_norm = relax_tol_correction_norm
 s% relax_tol_max_correction = relax_tol_max_correction
 s% relax_solver_iters_timestep_limit = relax_solver_iters_timestep_limit
 s% relax_iter_for_resid_tol2 = relax_iter_for_resid_tol2
 s% relax_tol_residual_norm1 = relax_tol_residual_norm1
 s% relax_tol_max_residual1 = relax_tol_max_residual1
 s% relax_iter_for_resid_tol3 = relax_iter_for_resid_tol3
 s% relax_tol_residual_norm2 = relax_tol_residual_norm2
 s% relax_tol_max_residual2 = relax_tol_max_residual2
 s% relax_tol_residual_norm3 = relax_tol_residual_norm3
 s% relax_tol_max_residual3 = relax_tol_max_residual3
 s% relax_maxT_for_gold_tolerances = relax_maxT_for_gold_tolerances
 
 s% use_gold_tolerances = use_gold_tolerances
 s% gold_solver_iters_timestep_limit = gold_solver_iters_timestep_limit 
 s% maxT_for_gold_tolerances = maxT_for_gold_tolerances
 s% gold_tol_residual_norm1 = gold_tol_residual_norm1
 s% gold_tol_max_residual1 = gold_tol_max_residual1
 s% gold_iter_for_resid_tol2 = gold_iter_for_resid_tol2
 s% gold_tol_residual_norm2 = gold_tol_residual_norm2
 s% gold_tol_max_residual2 = gold_tol_max_residual2
 s% gold_iter_for_resid_tol3 = gold_iter_for_resid_tol3
 s% gold_tol_residual_norm3 = gold_tol_residual_norm3
 s% gold_tol_max_residual3 = gold_tol_max_residual3
 s% steps_before_use_gold_tolerances = steps_before_use_gold_tolerances
 
 s% use_gold2_tolerances = use_gold2_tolerances
 s% gold2_solver_iters_timestep_limit = gold2_solver_iters_timestep_limit 
 s% gold2_tol_residual_norm1 = gold2_tol_residual_norm1
 s% gold2_tol_max_residual1 = gold2_tol_max_residual1
 s% gold2_iter_for_resid_tol2 = gold2_iter_for_resid_tol2
 s% gold2_tol_residual_norm2 = gold2_tol_residual_norm2
 s% gold2_tol_max_residual2 = gold2_tol_max_residual2
 s% gold2_iter_for_resid_tol3 = gold2_iter_for_resid_tol3
 s% gold2_tol_residual_norm3 = gold2_tol_residual_norm3
 s% gold2_tol_max_residual3 = gold2_tol_max_residual3
 s% steps_before_use_gold2_tolerances = steps_before_use_gold2_tolerances
 
 s% include_rotation_in_total_energy = include_rotation_in_total_energy

 s% convergence_ignore_equL_residuals = convergence_ignore_equL_residuals
 s% convergence_ignore_alpha_RTI_residuals = convergence_ignore_alpha_RTI_residuals

 s% iter_for_resid_tol2 = iter_for_resid_tol2
 s% iter_for_resid_tol3 = iter_for_resid_tol3

 s% solver_itermin = solver_itermin
 s% solver_itermin_until_reduce_min_corr_coeff = solver_itermin_until_reduce_min_corr_coeff
 s% solver_reduced_min_corr_coeff = solver_reduced_min_corr_coeff
 s% do_solver_damping_for_neg_xa = do_solver_damping_for_neg_xa
 s% hydro_mtx_max_allowed_abs_dlogT = hydro_mtx_max_allowed_abs_dlogT
 s% hydro_mtx_max_allowed_abs_dlogRho = hydro_mtx_max_allowed_abs_dlogRho
 s% min_logT_for_hydro_mtx_max_allowed = min_logT_for_hydro_mtx_max_allowed
 s% hydro_mtx_max_allowed_logT = hydro_mtx_max_allowed_logT
 s% hydro_mtx_max_allowed_logRho = hydro_mtx_max_allowed_logRho
 s% hydro_mtx_min_allowed_logT = hydro_mtx_min_allowed_logT
 s% hydro_mtx_min_allowed_logRho = hydro_mtx_min_allowed_logRho
 
 s% use_DGESVX_in_bcyclic = use_DGESVX_in_bcyclic
 s% use_equilibration_in_DGESVX = use_equilibration_in_DGESVX
 s% report_min_rcond_from_DGESXV = report_min_rcond_from_DGESXV
 
 s% op_split_burn = op_split_burn
 s% op_split_burn_min_T = op_split_burn_min_T
 s% op_split_burn_eps = op_split_burn_eps
 s% op_split_burn_odescal = op_split_burn_odescal
 s% op_split_burn_min_T_for_variable_T_solver = op_split_burn_min_T_for_variable_T_solver

 s% tiny_corr_coeff_limit = tiny_corr_coeff_limit
 s% scale_correction_norm = scale_correction_norm
 s% num_times_solver_reuse_mtx = num_times_solver_reuse_mtx
 s% corr_param_factor = corr_param_factor
 s% scale_max_correction = scale_max_correction
 s% ignore_min_corr_coeff_for_scale_max_correction = ignore_min_corr_coeff_for_scale_max_correction
 s% ignore_too_large_correction = ignore_too_large_correction
 s% ignore_species_in_max_correction = ignore_species_in_max_correction

 s% corr_norm_jump_limit = corr_norm_jump_limit
 s% max_corr_jump_limit = max_corr_jump_limit
 s% resid_norm_jump_limit = resid_norm_jump_limit
 s% max_resid_jump_limit = max_resid_jump_limit

 s% corr_coeff_limit = corr_coeff_limit
 s% tiny_corr_factor = tiny_corr_factor

 s% solver_max_tries_before_reject = solver_max_tries_before_reject
 s% max_tries1 = max_tries1
 s% max_tries_for_retry = max_tries_for_retry
 s% max_tries_after_5_retries = max_tries_after_5_retries
 s% max_tries_after_10_retries = max_tries_after_10_retries
 s% max_tries_after_20_retries = max_tries_after_20_retries
 s% retry_limit = retry_limit
 s% redo_limit = redo_limit

 s% use_Pvsc_art_visc = use_Pvsc_art_visc
 s% Pvsc_cq = Pvsc_cq
 s% Pvsc_zsh = Pvsc_zsh

 s% min_xa_hard_limit = min_xa_hard_limit
 s% min_xa_hard_limit_for_highT = min_xa_hard_limit_for_highT
 s% logT_max_for_min_xa_hard_limit = logT_max_for_min_xa_hard_limit
 s% logT_min_for_min_xa_hard_limit_for_highT = logT_min_for_min_xa_hard_limit_for_highT

 s% sum_xa_hard_limit = sum_xa_hard_limit
 s% sum_xa_hard_limit_for_highT = sum_xa_hard_limit_for_highT
 s% logT_max_for_sum_xa_hard_limit = logT_max_for_sum_xa_hard_limit
 s% logT_min_for_sum_xa_hard_limit_for_highT = logT_min_for_sum_xa_hard_limit_for_highT

 s% xa_clip_limit = xa_clip_limit
 s% report_solver_progress = report_solver_progress
 s% solver_test_partials_call_number = solver_test_partials_call_number
 s% solver_test_partials_iter_number = solver_test_partials_iter_number
 s% solver_epsder_chem = solver_epsder_chem
 s% solver_epsder_struct = solver_epsder_struct
 s% solver_numerical_jacobian = solver_numerical_jacobian
 s% solver_jacobian_nzlo = solver_jacobian_nzlo
 s% solver_jacobian_nzhi = solver_jacobian_nzhi
 s% solver_check_everything = solver_check_everything
 s% energy_conservation_dump_model_number = energy_conservation_dump_model_number
 s% solver_inspect_soln_flag = solver_inspect_soln_flag
 s% solver_test_partials_dx_0 = solver_test_partials_dx_0
 s% solver_test_partials_k = solver_test_partials_k
 s% solver_test_partials_k_low = solver_test_partials_k_low
 s% solver_test_partials_k_high = solver_test_partials_k_high
 s% solver_show_correction_info = solver_show_correction_info
 s% solver_test_partials_write_eos_call_info = solver_test_partials_write_eos_call_info
 s% solver_test_eos_partials = solver_test_eos_partials
 s% solver_test_kap_partials = solver_test_kap_partials
 s% solver_test_net_partials = solver_test_net_partials
 s% solver_test_atm_partials = solver_test_atm_partials
 s% solver_test_partials_var_name = solver_test_partials_var_name
 s% solver_test_partials_sink_name = solver_test_partials_sink_name
 s% solver_test_partials_equ_name = solver_test_partials_equ_name
 s% solver_test_partials_show_dx_var_name = solver_test_partials_show_dx_var_name
 s% solver_save_photo_call_number = solver_save_photo_call_number
 s% fill_arrays_with_NaNs = fill_arrays_with_NaNs
 s% zero_when_allocate = zero_when_allocate
 s% warn_when_large_rel_run_E_err = warn_when_large_rel_run_E_err
 s% warn_when_large_virial_thm_rel_err = warn_when_large_virial_thm_rel_err
 s% warn_when_get_a_bad_eos_result = warn_when_get_a_bad_eos_result
 s% warn_rates_for_high_temp = warn_rates_for_high_temp
 s% max_safe_logT_for_rates = max_safe_logT_for_rates
 s% eps_mdot_leak_frac_factor = eps_mdot_leak_frac_factor

 s% max_dt_div_tau_conv_for_TDC = max_dt_div_tau_conv_for_TDC
 s% max_dt_years_for_TDC = max_dt_years_for_TDC
 s% alpha_TDC_DAMP = alpha_TDC_DAMP
 s% alpha_TDC_DAMPR = alpha_TDC_DAMPR
 s% alpha_TDC_PtdVdt = alpha_TDC_PtdVdt
 s% max_X_for_gradT_eqn = max_X_for_gradT_eqn
 s% compare_TDC_to_MLT = compare_TDC_to_MLT

 s% RSP2_alfap = RSP2_alfap
 s% RSP2_alfad = RSP2_alfad
 s% RSP2_alfat = RSP2_alfat 
 s% RSP2_alfam = RSP2_alfam
 s% RSP2_alfar = RSP2_alfar
 s% RSP2_min_Lt_div_L_for_overshooting_mixing_type = RSP2_min_Lt_div_L_for_overshooting_mixing_type
 s% RSP2_min_Lc_div_L_for_convective_mixing_type = RSP2_min_Lc_div_L_for_convective_mixing_type
 s% RSP2_Lsurf_factor = RSP2_Lsurf_factor
 s% RSP2_use_Stellingwerf_Lr = RSP2_use_Stellingwerf_Lr
 s% RSP2_use_L_eqn_at_surface = RSP2_use_L_eqn_at_surface
 s% RSP2_assume_HSE = RSP2_assume_HSE
 s% RSP2_use_RSP_eqn_for_Y_face = RSP2_use_RSP_eqn_for_Y_face
 s% RSP2_use_mass_interp_face_values = RSP2_use_mass_interp_face_values
 s% RSP2_num_outermost_cells_forced_nonturbulent = RSP2_num_outermost_cells_forced_nonturbulent
 s% RSP2_num_innermost_cells_forced_nonturbulent = RSP2_num_innermost_cells_forced_nonturbulent
 s% RSP2_min_dt_div_tau_conv_switch_to_MLT = RSP2_min_dt_div_tau_conv_switch_to_MLT
 s% RSP2_min_dt_years_switch_to_MLT = RSP2_min_dt_years_switch_to_MLT
 s% RSP2_target_steps_per_cycle = RSP2_target_steps_per_cycle
 s% RSP2_max_num_periods = RSP2_max_num_periods
 s% RSP2_work_period = RSP2_work_period
 s% RSP2_map_first_period = RSP2_map_first_period
 s% RSP2_map_last_period = RSP2_map_last_period
 s% RSP2_min_max_R_for_periods = RSP2_min_max_R_for_periods
 s% RSP2_GREKM_avg_abs_frac_new = RSP2_GREKM_avg_abs_frac_new
 s% RSP2_GREKM_avg_abs_limit = RSP2_GREKM_avg_abs_limit
 s% RSP2_map_zone_interval = RSP2_map_zone_interval
 s% RSP2_work_filename = RSP2_work_filename
 s% RSP2_map_columns_filename = RSP2_map_columns_filename
 s% RSP2_map_filename = RSP2_map_filename
 s% RSP2_map_history_filename = RSP2_map_history_filename
 s% RSP2_write_map = RSP2_write_map
 s% RSP2_w_min_for_damping = RSP2_w_min_for_damping
 s% RSP2_source_seed = RSP2_source_seed
 s% RSP2_w_fix_if_neg = RSP2_w_fix_if_neg
 
 s% max_X_for_conv_timescale = max_X_for_conv_timescale
 s% min_X_for_conv_timescale = min_X_for_conv_timescale
 s% max_q_for_conv_timescale = max_q_for_conv_timescale
 s% min_q_for_conv_timescale = min_q_for_conv_timescale
 s% max_q_for_QHSE_timescale = max_q_for_QHSE_timescale
 s% min_q_for_QHSE_timescale = min_q_for_QHSE_timescale

 ! timestep
 s% max_timestep = max_timestep
 s% max_years_for_timestep = max_years_for_timestep

 s% hi_T_max_years_for_timestep = hi_T_max_years_for_timestep
 s% max_timestep_hi_T_limit = max_timestep_hi_T_limit

 s% min_timestep_factor = min_timestep_factor
 s% max_timestep_factor = max_timestep_factor
 s% max_timestep_factor_at_high_T = max_timestep_factor_at_high_T
 s% min_logT_for_max_timestep_factor_at_high_T = min_logT_for_max_timestep_factor_at_high_T
 s% time_delta_coeff = time_delta_coeff
 s% timestep_factor_for_retries = timestep_factor_for_retries
 s% retry_hold = retry_hold
 s% neg_mass_fraction_hold = neg_mass_fraction_hold
 s% timestep_dt_factor = timestep_dt_factor
 s% use_dt_low_pass_controller = use_dt_low_pass_controller
 
 s% force_timestep_min = force_timestep_min
 s% force_timestep_min_years = force_timestep_min_years
 s% force_timestep_min_factor = force_timestep_min_factor
 s% force_timestep = force_timestep
 s% force_timestep_years = force_timestep_years

 s% varcontrol_target = varcontrol_target
 s% min_allowed_varcontrol_target = min_allowed_varcontrol_target
 s% varcontrol_dt_limit_ratio_hard_max = varcontrol_dt_limit_ratio_hard_max
 s% xa_scale = xa_scale

 s% solver_iters_timestep_limit = solver_iters_timestep_limit

 s% burn_steps_limit = burn_steps_limit
 s% burn_steps_hard_limit = burn_steps_hard_limit

 s% diffusion_steps_limit = diffusion_steps_limit
 s% diffusion_steps_hard_limit = diffusion_steps_hard_limit
 s% diffusion_iters_limit = diffusion_iters_limit
 s% diffusion_iters_hard_limit = diffusion_iters_hard_limit

 s% dt_div_dt_cell_collapse_limit = dt_div_dt_cell_collapse_limit
 s% dt_div_dt_cell_collapse_hard_limit = dt_div_dt_cell_collapse_hard_limit
 s% dt_div_min_dr_div_cs_limit = dt_div_min_dr_div_cs_limit
 s% dt_div_min_dr_div_cs_hard_limit = dt_div_min_dr_div_cs_hard_limit
 
 s% min_abs_du_div_cs_for_dt_div_min_dr_div_cs_limit = min_abs_du_div_cs_for_dt_div_min_dr_div_cs_limit
 s% min_abs_u_div_cs_for_dt_div_min_dr_div_cs_limit = min_abs_u_div_cs_for_dt_div_min_dr_div_cs_limit
 s% min_k_for_dt_div_min_dr_div_cs_limit = min_k_for_dt_div_min_dr_div_cs_limit
 s% min_q_for_dt_div_min_dr_div_cs_limit = min_q_for_dt_div_min_dr_div_cs_limit
 s% max_q_for_dt_div_min_dr_div_cs_limit = max_q_for_dt_div_min_dr_div_cs_limit
 s% check_remnant_only_for_dt_div_min_dr_div_cs_limit = check_remnant_only_for_dt_div_min_dr_div_cs_limit

 s% dX_mix_dist_limit = dX_mix_dist_limit

 s% dH_limit_min_H = dH_limit_min_H
 s% dH_limit = dH_limit
 s% dH_hard_limit = dH_hard_limit
 s% dH_div_H_limit_min_H = dH_div_H_limit_min_H
 s% dH_div_H_limit = dH_div_H_limit
 s% dH_div_H_hard_limit = dH_div_H_hard_limit
 s% dH_decreases_only = dH_decreases_only

 s% dHe_limit_min_He = dHe_limit_min_He
 s% dHe_limit = dHe_limit
 s% dHe_hard_limit = dHe_hard_limit
 s% dHe_div_He_limit_min_He = dHe_div_He_limit_min_He
 s% dHe_div_He_limit = dHe_div_He_limit
 s% dHe_div_He_hard_limit = dHe_div_He_hard_limit
 s% dHe_decreases_only = dHe_decreases_only

 s% dHe3_limit_min_He3 = dHe3_limit_min_He3
 s% dHe3_limit = dHe3_limit
 s% dHe3_hard_limit = dHe3_hard_limit
 s% dHe3_div_He3_limit_min_He3 = dHe3_div_He3_limit_min_He3
 s% dHe3_div_He3_limit = dHe3_div_He3_limit
 s% dHe3_div_He3_hard_limit = dHe3_div_He3_hard_limit
 s% dHe3_decreases_only = dHe3_decreases_only

 s% dX_limit_min_X = dX_limit_min_X
 s% dX_limit = dX_limit
 s% dX_hard_limit = dX_hard_limit
 s% dX_div_X_limit_min_X = dX_div_X_limit_min_X
 s% dX_div_X_limit = dX_div_X_limit
 s% dX_div_X_hard_limit = dX_div_X_hard_limit
 s% dX_div_X_at_high_T_limit = dX_div_X_at_high_T_limit
 s% dX_div_X_at_high_T_hard_limit = dX_div_X_at_high_T_hard_limit
 s% dX_div_X_at_high_T_limit_lgT_min = dX_div_X_at_high_T_limit_lgT_min
 
 s% dX_decreases_only = dX_decreases_only

 s% dX_nuc_drop_min_X_limit = dX_nuc_drop_min_X_limit
 s% dX_nuc_drop_max_A_limit = dX_nuc_drop_max_A_limit
 s% dX_nuc_drop_limit = dX_nuc_drop_limit
 s% dX_nuc_drop_limit_at_high_T = dX_nuc_drop_limit_at_high_T
 s% dX_nuc_drop_hard_limit = dX_nuc_drop_hard_limit
 s% dX_nuc_drop_min_yrs_for_dt = dX_nuc_drop_min_yrs_for_dt

 s% dL_div_L_limit_min_L = dL_div_L_limit_min_L
 s% dL_div_L_limit = dL_div_L_limit
 s% dL_div_L_hard_limit = dL_div_L_hard_limit

 s% delta_lgP_limit = delta_lgP_limit
 s% delta_lgP_hard_limit = delta_lgP_hard_limit
 s% delta_lgP_limit_min_lgP = delta_lgP_limit_min_lgP

 s% delta_lgRho_limit = delta_lgRho_limit
 s% delta_lgRho_hard_limit = delta_lgRho_hard_limit
 s% delta_lgRho_limit_min_lgRho = delta_lgRho_limit_min_lgRho

 s% delta_lgT_limit = delta_lgT_limit
 s% delta_lgT_hard_limit = delta_lgT_hard_limit
 s% delta_lgT_limit_min_lgT = delta_lgT_limit_min_lgT

 s% delta_lgE_limit = delta_lgE_limit
 s% delta_lgE_hard_limit = delta_lgE_hard_limit
 s% delta_lgE_limit_min_lgE = delta_lgE_limit_min_lgE

 s% delta_lgR_limit = delta_lgR_limit
 s% delta_lgR_hard_limit = delta_lgR_hard_limit
 s% delta_lgR_limit_min_lgR = delta_lgR_limit_min_lgR

 s% delta_Ye_highT_limit = delta_Ye_highT_limit
 s% delta_Ye_highT_hard_limit = delta_Ye_highT_hard_limit
 s% minT_for_highT_Ye_limit = minT_for_highT_Ye_limit

 s% delta_lgL_nuc_cat_limit = delta_lgL_nuc_cat_limit
 s% delta_lgL_nuc_cat_hard_limit = delta_lgL_nuc_cat_hard_limit
 s% lgL_nuc_cat_burn_min = lgL_nuc_cat_burn_min
 s% lgL_nuc_mix_dist_limit = lgL_nuc_mix_dist_limit

 s% delta_lgL_H_limit = delta_lgL_H_limit
 s% delta_lgL_H_hard_limit = delta_lgL_H_hard_limit
 s% lgL_H_burn_min = lgL_H_burn_min
 s% lgL_H_drop_factor = lgL_H_drop_factor
 s% lgL_H_burn_relative_limit = lgL_H_burn_relative_limit

 s% delta_lgL_He_limit = delta_lgL_He_limit
 s% delta_lgL_He_hard_limit = delta_lgL_He_hard_limit
 s% lgL_He_burn_min = lgL_He_burn_min
 s% lgL_He_drop_factor = lgL_He_drop_factor
 s% lgL_He_burn_relative_limit = lgL_He_burn_relative_limit

 s% delta_lgL_z_limit = delta_lgL_z_limit
 s% delta_lgL_z_hard_limit = delta_lgL_z_hard_limit
 s% lgL_z_burn_min = lgL_z_burn_min
 s% lgL_z_drop_factor = lgL_z_drop_factor
 s% lgL_z_burn_relative_limit = lgL_z_burn_relative_limit

 s% delta_lgL_power_photo_limit = delta_lgL_power_photo_limit
 s% delta_lgL_power_photo_hard_limit = delta_lgL_power_photo_hard_limit
 s% lgL_power_photo_burn_min = lgL_power_photo_burn_min
 s% lgL_power_photo_drop_factor = lgL_power_photo_drop_factor
 s% min_lgT_for_lgL_power_photo_limit = min_lgT_for_lgL_power_photo_limit

 s% delta_lgL_nuc_limit = delta_lgL_nuc_limit
 s% delta_lgL_nuc_hard_limit = delta_lgL_nuc_hard_limit
 s% delta_lgL_nuc_at_high_T_limit = delta_lgL_nuc_at_high_T_limit
 s% delta_lgL_nuc_at_high_T_hard_limit = delta_lgL_nuc_at_high_T_hard_limit
 s% delta_lgL_nuc_at_high_T_limit_lgT_min = delta_lgL_nuc_at_high_T_limit_lgT_min
 
 s% max_lgT_for_lgL_nuc_limit = max_lgT_for_lgL_nuc_limit
 s% lgL_nuc_burn_min = lgL_nuc_burn_min
 s% lgL_nuc_drop_factor = lgL_nuc_drop_factor

 s% delta_lgRho_cntr_limit = delta_lgRho_cntr_limit
 s% delta_lgRho_cntr_hard_limit = delta_lgRho_cntr_hard_limit

 s% delta_lgP_cntr_limit = delta_lgP_cntr_limit
 s% delta_lgP_cntr_hard_limit = delta_lgP_cntr_hard_limit

 s% delta_lgT_cntr_limit = delta_lgT_cntr_limit
 s% delta_lgT_cntr_hard_limit = delta_lgT_cntr_hard_limit
 s% delta_lgT_cntr_limit_only_after_near_zams = delta_lgT_cntr_limit_only_after_near_zams

 s% delta_lgT_max_limit = delta_lgT_max_limit
 s% delta_lgT_max_hard_limit = delta_lgT_max_hard_limit
 s% delta_lgT_max_limit_lgT_min = delta_lgT_max_limit_lgT_min
 s% delta_lgT_max_limit_only_after_near_zams = delta_lgT_max_limit_only_after_near_zams

 s% delta_lgT_max_at_high_T_limit = delta_lgT_max_at_high_T_limit
 s% delta_lgT_max_at_high_T_hard_limit = delta_lgT_max_at_high_T_hard_limit
 s% delta_lgT_max_at_high_T_limit_lgT_min = delta_lgT_max_at_high_T_limit_lgT_min

 s% delta_log_eps_nuc_limit = delta_log_eps_nuc_limit
 s% delta_log_eps_nuc_hard_limit = delta_log_eps_nuc_hard_limit

 s% delta_dX_div_X_cntr_min = delta_dX_div_X_cntr_min
 s% delta_dX_div_X_cntr_max = delta_dX_div_X_cntr_max
 s% delta_dX_div_X_cntr_limit = delta_dX_div_X_cntr_limit
 s% delta_dX_div_X_cntr_hard_limit = delta_dX_div_X_cntr_hard_limit

 s% delta_dX_div_X_drop_only = delta_dX_div_X_drop_only
 s% delta_lg_XH_drop_only = delta_lg_XH_drop_only
 s% delta_lg_XHe_drop_only = delta_lg_XHe_drop_only
 s% delta_lg_XC_drop_only = delta_lg_XC_drop_only
 s% delta_lg_XNe_drop_only = delta_lg_XNe_drop_only
 s% delta_lg_XO_drop_only = delta_lg_XO_drop_only
 s% delta_lg_XSi_drop_only = delta_lg_XSi_drop_only
 s% delta_XH_drop_only = delta_XH_drop_only
 s% delta_XHe_drop_only = delta_XHe_drop_only
 s% delta_XC_drop_only = delta_XC_drop_only
 s% delta_XNe_drop_only = delta_XNe_drop_only
 s% delta_XO_drop_only = delta_XO_drop_only
 s% delta_XSi_drop_only = delta_XSi_drop_only

 s% delta_lg_XH_cntr_min = delta_lg_XH_cntr_min
 s% delta_lg_XH_cntr_max = delta_lg_XH_cntr_max
 s% delta_lg_XH_cntr_limit = delta_lg_XH_cntr_limit
 s% delta_lg_XH_cntr_hard_limit = delta_lg_XH_cntr_hard_limit

 s% delta_lg_XHe_cntr_min = delta_lg_XHe_cntr_min
 s% delta_lg_XHe_cntr_max = delta_lg_XHe_cntr_max
 s% delta_lg_XHe_cntr_limit = delta_lg_XHe_cntr_limit
 s% delta_lg_XHe_cntr_hard_limit = delta_lg_XHe_cntr_hard_limit

 s% delta_lg_XC_cntr_min = delta_lg_XC_cntr_min
 s% delta_lg_XC_cntr_max = delta_lg_XC_cntr_max
 s% delta_lg_XC_cntr_limit = delta_lg_XC_cntr_limit
 s% delta_lg_XC_cntr_hard_limit = delta_lg_XC_cntr_hard_limit

 s% delta_lg_XNe_cntr_limit = delta_lg_XNe_cntr_limit
 s% delta_lg_XNe_cntr_hard_limit = delta_lg_XNe_cntr_hard_limit
 s% delta_lg_XNe_cntr_min = delta_lg_XNe_cntr_min
 s% delta_lg_XNe_cntr_max = delta_lg_XNe_cntr_max

 s% delta_lg_XO_cntr_limit = delta_lg_XO_cntr_limit
 s% delta_lg_XO_cntr_hard_limit = delta_lg_XO_cntr_hard_limit
 s% delta_lg_XO_cntr_min = delta_lg_XO_cntr_min
 s% delta_lg_XO_cntr_max = delta_lg_XO_cntr_max

 s% delta_lg_XSi_cntr_limit = delta_lg_XSi_cntr_limit
 s% delta_lg_XSi_cntr_hard_limit = delta_lg_XSi_cntr_hard_limit
 s% delta_lg_XSi_cntr_min = delta_lg_XSi_cntr_min
 s% delta_lg_XSi_cntr_max = delta_lg_XSi_cntr_max

 s% delta_XH_cntr_limit = delta_XH_cntr_limit
 s% delta_XH_cntr_hard_limit = delta_XH_cntr_hard_limit
 s% delta_XHe_cntr_limit = delta_XHe_cntr_limit
 s% delta_XHe_cntr_hard_limit = delta_XHe_cntr_hard_limit
 s% delta_XC_cntr_limit = delta_XC_cntr_limit
 s% delta_XC_cntr_hard_limit = delta_XC_cntr_hard_limit
 s% delta_XNe_cntr_limit = delta_XNe_cntr_limit
 s% delta_XNe_cntr_hard_limit = delta_XNe_cntr_hard_limit
 s% delta_XO_cntr_limit = delta_XO_cntr_limit
 s% delta_XO_cntr_hard_limit = delta_XO_cntr_hard_limit
 s% delta_XSi_cntr_limit = delta_XSi_cntr_limit
 s% delta_XSi_cntr_hard_limit = delta_XSi_cntr_hard_limit

 s% delta_lgTeff_limit = delta_lgTeff_limit
 s% delta_lgTeff_hard_limit = delta_lgTeff_hard_limit

 s% delta_lgL_limit = delta_lgL_limit
 s% delta_lgL_limit_L_min = delta_lgL_limit_L_min
 s% delta_lgL_hard_limit = delta_lgL_hard_limit

 s% delta_HR_ds_L = delta_HR_ds_L
 s% delta_HR_ds_Teff = delta_HR_ds_Teff
 s% delta_HR_limit = delta_HR_limit
 s% delta_HR_hard_limit = delta_HR_hard_limit

 s% delta_lg_star_mass_limit = delta_lg_star_mass_limit
 s% delta_lg_star_mass_hard_limit = delta_lg_star_mass_hard_limit

 s% delta_mdot_atol = delta_mdot_atol
 s% delta_mdot_rtol = delta_mdot_rtol
 s% delta_mdot_limit = delta_mdot_limit
 s% delta_mdot_hard_limit = delta_mdot_hard_limit

 s% adjust_J_q_limit = adjust_J_q_limit
 s% adjust_J_q_hard_limit = adjust_J_q_hard_limit
 s% never_skip_hard_limits = never_skip_hard_limits
 s% relax_hard_limits_after_retry = relax_hard_limits_after_retry
 s% report_dt_hard_limit_retries = report_dt_hard_limit_retries
 s% report_min_dr_div_cs = report_min_dr_div_cs
 s% report_solver_dt_info = report_solver_dt_info

 s% limit_for_rel_error_in_energy_conservation = limit_for_rel_error_in_energy_conservation
 s% hard_limit_for_rel_error_in_energy_conservation = hard_limit_for_rel_error_in_energy_conservation

 s% min_chem_eqn_scale = min_chem_eqn_scale

 s% trace_evolve = trace_evolve


 ! misc
 s% zams_filename = zams_filename
 s% set_rho_to_dm_div_dV = set_rho_to_dm_div_dV

 s% use_other_eos = use_other_eos
 s% use_other_surface_PT = use_other_surface_PT
 s% use_other_kap = use_other_kap
 s% use_other_diffusion = use_other_diffusion
 s% use_other_diffusion_factor = use_other_diffusion_factor
 s% use_other_adjust_mdot = use_other_adjust_mdot
 s% use_other_j_for_adjust_J_lost = use_other_j_for_adjust_J_lost
 s% use_other_alpha_mlt = use_other_alpha_mlt
 s% use_other_am_mixing = use_other_am_mixing
 s% use_other_brunt = use_other_brunt
 s% use_other_brunt_smoothing = use_other_brunt_smoothing
 s% use_other_solver_monitor = use_other_solver_monitor
 s% use_other_build_initial_model = use_other_build_initial_model
 s% use_other_cgrav = use_other_cgrav
 s% use_other_mesh_delta_coeff_factor = use_other_mesh_delta_coeff_factor
 s% use_other_energy_implicit = use_other_energy_implicit
 s% use_other_remove_surface = use_other_remove_surface
 s% use_other_momentum = use_other_momentum
 s% use_other_momentum_implicit = use_other_momentum_implicit
 s% use_other_pressure = use_other_pressure
 s% use_other_energy = use_other_energy
 s% use_other_mesh_functions = use_other_mesh_functions
 s% use_other_eps_grav = use_other_eps_grav
 s% use_other_gradr_factor = use_other_gradr_factor
 s% use_other_D_mix = use_other_D_mix
 s% use_other_neu = use_other_neu
 s% use_other_net_get = use_other_net_get
 s% use_other_opacity_factor = use_other_opacity_factor
 s% use_other_diffusion_coefficients = use_other_diffusion_coefficients
 s% use_other_pgstar_plots = use_other_pgstar_plots
 s% use_other_eval_fp_ft = use_other_eval_fp_ft
 s% use_other_eval_i_rot = use_other_eval_i_rot
 s% use_other_torque = use_other_torque
 s% use_other_torque_implicit = use_other_torque_implicit
 s% use_other_wind = use_other_wind
 s% use_other_accreting_state = use_other_accreting_state
 s% use_other_after_struct_burn_mix = use_other_after_struct_burn_mix
 s% use_other_before_struct_burn_mix = use_other_before_struct_burn_mix
 s% use_other_astero_freq_corr = use_other_astero_freq_corr
 s% use_other_timestep_limit = use_other_timestep_limit
 s% use_other_set_pgstar_controls = use_other_set_pgstar_controls
 s% use_other_screening = use_other_screening

 s% x_ctrl = x_ctrl
 s% x_integer_ctrl = x_integer_ctrl
 s% x_logical_ctrl = x_logical_ctrl
 s% x_character_ctrl = x_character_ctrl

 ! info for debugging
 s% stop_for_bad_nums = stop_for_bad_nums
 s% report_ierr = report_ierr
 s% report_bad_negative_xa = report_bad_negative_xa

 s% diffusion_dump_call_number = diffusion_dump_call_number


 end subroutine store_controls


 subroutine set_controls_for_writing(s, ierr)
 use star_private_def
 use chem_def ! categories
 type (star_info), pointer :: s
 integer, intent(out) :: ierr

 ierr = 0

 ! where to start
 initial_mass = s% initial_mass
 initial_z = s% initial_z
 initial_y = s% initial_y
 initial_he3 = s% initial_he3

 ! definition of core boundaries
 he_core_boundary_h1_fraction = s% he_core_boundary_h1_fraction
 co_core_boundary_he4_fraction = s% co_core_boundary_he4_fraction
 one_core_boundary_he4_c12_fraction = s% one_core_boundary_he4_c12_fraction
 fe_core_boundary_si28_fraction = s% fe_core_boundary_si28_fraction
 neutron_rich_core_boundary_Ye_max = s% neutron_rich_core_boundary_Ye_max
 min_boundary_fraction = s% min_boundary_fraction

 ! when to stop
 max_model_number = s% max_model_number
 max_abs_rel_run_E_err = s% max_abs_rel_run_E_err
 max_number_retries = s% max_number_retries
 relax_max_number_retries = s% relax_max_number_retries
 max_age = s% max_age
 max_age_in_days = s% max_age_in_days
 max_age_in_seconds = s% max_age_in_seconds
 num_adjusted_dt_steps_before_max_age = s% num_adjusted_dt_steps_before_max_age
 dt_years_for_steps_before_max_age = s% dt_years_for_steps_before_max_age
 reduction_factor_for_max_timestep = s% reduction_factor_for_max_timestep
 when_to_stop_rtol = s% when_to_stop_rtol
 when_to_stop_atol = s% when_to_stop_atol
 gamma_center_limit = s% gamma_center_limit
 eta_center_limit = s% eta_center_limit
 log_center_temp_limit = s% log_center_temp_limit
 log_max_temp_upper_limit = s% log_max_temp_upper_limit
 log_max_temp_lower_limit = s% log_max_temp_lower_limit
 log_center_temp_lower_limit = s% log_center_temp_lower_limit
 log_center_density_limit = s% log_center_density_limit
 log_center_density_lower_limit = s% log_center_density_lower_limit
 min_timestep_limit = s% min_timestep_limit

 center_entropy_limit = s% center_entropy_limit
 center_entropy_lower_limit = s% center_entropy_lower_limit
 max_entropy_limit = s% max_entropy_limit
 max_entropy_lower_limit = s% max_entropy_lower_limit

 fe_core_infall_limit = s% fe_core_infall_limit
 center_Ye_lower_limit = s% center_Ye_lower_limit
 center_R_lower_limit = s% center_R_lower_limit
 non_fe_core_infall_limit = s% non_fe_core_infall_limit
 non_fe_core_rebound_limit = s% non_fe_core_rebound_limit
 v_div_csound_surf_limit = s% v_div_csound_surf_limit
 v_div_csound_max_limit = s% v_div_csound_max_limit
 Lnuc_div_L_upper_limit = s% Lnuc_div_L_upper_limit
 Lnuc_div_L_lower_limit = s% Lnuc_div_L_lower_limit
 v_surf_div_v_kh_upper_limit = s% v_surf_div_v_kh_upper_limit
 v_surf_div_v_kh_lower_limit = s% v_surf_div_v_kh_lower_limit
 v_surf_div_v_esc_limit = s% v_surf_div_v_esc_limit
 v_surf_kms_limit = s% v_surf_kms_limit

 stop_near_zams = s% stop_near_zams
 stop_at_phase_PreMS = s% stop_at_phase_PreMS
 stop_at_phase_ZAMS = s% stop_at_phase_ZAMS
 stop_at_phase_IAMS = s% stop_at_phase_IAMS
 stop_at_phase_TAMS = s% stop_at_phase_TAMS
 stop_at_phase_He_Burn = s% stop_at_phase_He_Burn
 stop_at_phase_ZACHeB = s% stop_at_phase_ZACHeB
 stop_at_phase_TACHeB = s% stop_at_phase_TACHeB
 stop_at_phase_TP_AGB = s% stop_at_phase_TP_AGB
 stop_at_phase_C_Burn = s% stop_at_phase_C_Burn
 stop_at_phase_Ne_Burn = s% stop_at_phase_Ne_Burn
 stop_at_phase_O_Burn = s% stop_at_phase_O_Burn
 stop_at_phase_Si_Burn = s% stop_at_phase_Si_Burn
 stop_at_phase_WDCS = s% stop_at_phase_WDCS
 Lnuc_div_L_zams_limit = s% Lnuc_div_L_zams_limit
 Pgas_div_P_limit = s% Pgas_div_P_limit
 Pgas_div_P_limit_max_q = s% Pgas_div_P_limit_max_q
 gamma1_limit = s% gamma1_limit
 gamma1_limit_max_q = s% gamma1_limit_max_q
 gamma1_limit_max_v_div_vesc = s% gamma1_limit_max_v_div_vesc
 peak_burn_vconv_div_cs_limit = s% peak_burn_vconv_div_cs_limit
 omega_div_omega_crit_limit = s% omega_div_omega_crit_limit
 delta_nu_lower_limit = s% delta_nu_lower_limit
 delta_nu_upper_limit = s% delta_nu_upper_limit
 delta_Pg_lower_limit = s% delta_Pg_lower_limit
 delta_Pg_upper_limit = s% delta_Pg_upper_limit
 shock_mass_upper_limit = s% shock_mass_upper_limit
 mach1_mass_upper_limit = s% mach1_mass_upper_limit
 stop_when_reach_this_cumulative_extra_heating = s% stop_when_reach_this_cumulative_extra_heating

 xa_central_lower_limit_species = s% xa_central_lower_limit_species
 xa_central_lower_limit = s% xa_central_lower_limit

 xa_central_upper_limit_species = s% xa_central_upper_limit_species
 xa_central_upper_limit = s% xa_central_upper_limit

 xa_surface_lower_limit_species = s% xa_surface_lower_limit_species
 xa_surface_lower_limit = s% xa_surface_lower_limit

 xa_surface_upper_limit_species = s% xa_surface_upper_limit_species
 xa_surface_upper_limit = s% xa_surface_upper_limit

 xa_average_lower_limit_species = s% xa_average_lower_limit_species
 xa_average_lower_limit = s% xa_average_lower_limit

 xa_average_upper_limit_species = s% xa_average_upper_limit_species
 xa_average_upper_limit = s% xa_average_upper_limit

 HB_limit = s% HB_limit

 star_mass_max_limit = s% star_mass_max_limit
 star_mass_min_limit = s% star_mass_min_limit
 ejecta_mass_max_limit = s% ejecta_mass_max_limit
 remnant_mass_min_limit = s% remnant_mass_min_limit
 
 star_species_mass_min_limit = s% star_species_mass_min_limit
 star_species_mass_min_limit_iso = s% star_species_mass_min_limit_iso
 star_species_mass_max_limit = s% star_species_mass_max_limit
 star_species_mass_max_limit_iso = s% star_species_mass_max_limit_iso

 xmstar_min_limit = s% xmstar_min_limit
 xmstar_max_limit = s% xmstar_max_limit
 envelope_mass_limit = s% envelope_mass_limit
 envelope_fraction_left_limit = s% envelope_fraction_left_limit

 he_core_mass_limit = s% he_core_mass_limit
 co_core_mass_limit = s% co_core_mass_limit
 one_core_mass_limit = s% one_core_mass_limit
 fe_core_mass_limit = s% fe_core_mass_limit
 neutron_rich_core_mass_limit = s% neutron_rich_core_mass_limit

 he_layer_mass_lower_limit = s% he_layer_mass_lower_limit
 abs_diff_lg_LH_lg_Ls_limit = s% abs_diff_lg_LH_lg_Ls_limit
 Teff_upper_limit = s% Teff_upper_limit
 Teff_lower_limit = s% Teff_lower_limit
 photosphere_m_upper_limit = s% photosphere_m_upper_limit
 photosphere_m_lower_limit = s% photosphere_m_lower_limit
 photosphere_m_sub_M_center_limit = s% photosphere_m_sub_M_center_limit
 photosphere_r_upper_limit = s% photosphere_r_upper_limit
 photosphere_r_lower_limit = s% photosphere_r_lower_limit
 log_Teff_upper_limit = s% log_Teff_upper_limit
 log_Teff_lower_limit = s% log_Teff_lower_limit
 log_Tsurf_upper_limit = s% log_Tsurf_upper_limit
 log_Tsurf_lower_limit = s% log_Tsurf_lower_limit
 log_Rsurf_upper_limit = s% log_Rsurf_upper_limit
 log_Rsurf_lower_limit = s% log_Rsurf_lower_limit
 log_Psurf_upper_limit = s% log_Psurf_upper_limit
 log_Psurf_lower_limit = s% log_Psurf_lower_limit
 log_Dsurf_upper_limit = s% log_Dsurf_upper_limit
 log_Dsurf_lower_limit = s% log_Dsurf_lower_limit
 log_L_upper_limit = s% log_L_upper_limit
 log_L_lower_limit = s% log_L_lower_limit
 log_g_upper_limit = s% log_g_upper_limit
 log_g_lower_limit = s% log_g_lower_limit

 power_nuc_burn_upper_limit = s% power_nuc_burn_upper_limit
 power_h_burn_upper_limit = s% power_h_burn_upper_limit
 power_he_burn_upper_limit = s% power_he_burn_upper_limit
 power_z_burn_upper_limit = s% power_z_burn_upper_limit
 power_nuc_burn_lower_limit = s% power_nuc_burn_lower_limit
 power_h_burn_lower_limit = s% power_h_burn_lower_limit
 power_he_burn_lower_limit = s% power_he_burn_lower_limit
 power_z_burn_lower_limit = s% power_z_burn_lower_limit


 ! output of "snapshots" for restarts
 photo_interval = s% photo_interval
 photo_digits = s% photo_digits
 photo_directory = s% photo_directory
 ! output of history and profiles.
 do_history_file = s% do_history_file
 history_interval = s% history_interval

 write_header_frequency = s% write_header_frequency
 terminal_interval = s% terminal_interval
 terminal_show_age_units = s% terminal_show_age_units
 terminal_show_timestep_units = s% terminal_show_timestep_units
 terminal_show_log_dt = s% terminal_show_log_dt
 terminal_show_log_age = s% terminal_show_log_age
 extra_terminal_output_file = s% extra_terminal_output_file
 num_trace_history_values = s% num_trace_history_values
 trace_history_value_name = s% trace_history_value_name

 log_directory = s% log_directory

 star_history_name = s% star_history_name
 star_history_header_name = s% star_history_header_name
 star_history_dbl_format = s% star_history_dbl_format
 star_history_int_format = s% star_history_int_format
 star_history_txt_format = s% star_history_txt_format

 profiles_index_name = s% profiles_index_name
 profile_data_prefix = s% profile_data_prefix
 profile_data_suffix = s% profile_data_suffix
 profile_data_header_suffix = s% profile_data_header_suffix
 profile_int_format = s% profile_int_format
 profile_txt_format = s% profile_txt_format
 profile_dbl_format = s% profile_dbl_format
 profile_header_include_sys_details = s% profile_header_include_sys_details
 write_profiles_flag = s% write_profiles_flag
 profile_interval = s% profile_interval
 priority_profile_interval = s% priority_profile_interval
 profile_model = s% profile_model
 max_num_profile_models = s% max_num_profile_models
 max_num_profile_zones = s% max_num_profile_zones

 write_controls_info_with_profile = s% write_controls_info_with_profile
 controls_data_prefix = s% controls_data_prefix
 controls_data_suffix = s% controls_data_suffix

 write_pulse_data_with_profile = s% write_pulse_data_with_profile
 pulse_data_format = s% pulse_data_format
 add_atmosphere_to_pulse_data = s% add_atmosphere_to_pulse_data
 add_center_point_to_pulse_data = s% add_center_point_to_pulse_data
 keep_surface_point_for_pulse_data = s% keep_surface_point_for_pulse_data
 add_double_points_to_pulse_data = s% add_double_points_to_pulse_data
 interpolate_rho_for_pulse_data = s% interpolate_rho_for_pulse_data
 threshold_grad_mu_for_double_point = s% threshold_grad_mu_for_double_point
 max_number_of_double_points = s% max_number_of_double_points

 fgong_header = s% fgong_header
 fgong_ivers = s% fgong_ivers
 
 max_num_gyre_points = s% max_num_gyre_points
 format_for_OSC_data = s% format_for_OSC_data
 fgong_zero_A_inside_r = s% fgong_zero_A_inside_r
 use_other_export_pulse_data = s% use_other_export_pulse_data
 use_other_get_pulse_data = s% use_other_get_pulse_data
 use_other_edit_pulse_data = s% use_other_edit_pulse_data

 write_model_with_profile = s% write_model_with_profile
 model_data_prefix = s% model_data_prefix
 model_data_suffix = s% model_data_suffix

 mixing_D_limit_for_log = s% mixing_D_limit_for_log
 trace_mass_location = s% trace_mass_location
 min_tau_for_max_abs_v_location = s% min_tau_for_max_abs_v_location
 min_q_for_inner_mach1_location = s% min_q_for_inner_mach1_location
 max_q_for_outer_mach1_location = s% max_q_for_outer_mach1_location
 
 mass_depth_for_L_surf = s% mass_depth_for_L_surf
 conv_core_gap_dq_limit = s% conv_core_gap_dq_limit

 ! burn zone eps definitions for use in logs and profiles
 burn_min1 = s% burn_min1
 burn_min2 = s% burn_min2

 max_conv_vel_div_csound_maxq = s% max_conv_vel_div_csound_maxq
 width_for_limit_conv_vel = s% width_for_limit_conv_vel
 max_q_for_limit_conv_vel = s% max_q_for_limit_conv_vel
 max_mass_in_gm_for_limit_conv_vel = s% max_mass_in_gm_for_limit_conv_vel
 max_r_in_cm_for_limit_conv_vel = s% max_r_in_cm_for_limit_conv_vel

 ! for reported average values
 surface_avg_abundance_dq = s% surface_avg_abundance_dq
 center_avg_value_dq = s% center_avg_value_dq

 ! mixing parameters
 min_convective_gap = s% min_convective_gap
 min_thermohaline_gap = s% min_thermohaline_gap
 min_semiconvection_gap = s% min_semiconvection_gap
 min_thermohaline_dropout = s% min_thermohaline_dropout
 max_dropout_gradL_sub_grada = s% max_dropout_gradL_sub_grada
 remove_embedded_semiconvection = s% remove_embedded_semiconvection
 recalc_mix_info_after_evolve = s% recalc_mix_info_after_evolve
 remove_mixing_glitches = s% remove_mixing_glitches
 okay_to_remove_mixing_singleton = s% okay_to_remove_mixing_singleton
 prune_bad_cz_min_Hp_height = s% prune_bad_cz_min_Hp_height
 prune_bad_cz_min_log_eps_nuc = s% prune_bad_cz_min_log_eps_nuc
 redo_conv_for_dr_lt_mixing_length = s% redo_conv_for_dr_lt_mixing_length

 alpha_semiconvection = s% alpha_semiconvection
 semiconvection_option = s% semiconvection_option
 use_Ledoux_criterion = s% use_Ledoux_criterion
 num_cells_for_smooth_gradL_composition_term = s% num_cells_for_smooth_gradL_composition_term
 threshold_for_smooth_gradL_composition_term = s% threshold_for_smooth_gradL_composition_term
 clip_D_limit = s% clip_D_limit
 fix_eps_grav_transition_to_grid = s% fix_eps_grav_transition_to_grid

 okay_to_reduce_gradT_excess = s% okay_to_reduce_gradT_excess
 gradT_excess_f1 = s% gradT_excess_f1
 gradT_excess_f2 = s% gradT_excess_f2
 gradT_excess_max_center_h1 = s% gradT_excess_max_center_h1
 gradT_excess_min_center_he4 = s% gradT_excess_min_center_he4
 gradT_excess_max_logT = s% gradT_excess_max_logT
 gradT_excess_min_log_tau_full_on = s% gradT_excess_min_log_tau_full_on
 gradT_excess_max_log_tau_full_off = s% gradT_excess_max_log_tau_full_off
 gradT_excess_lambda1 = s% gradT_excess_lambda1
 gradT_excess_beta1 = s% gradT_excess_beta1
 gradT_excess_lambda2 = s% gradT_excess_lambda2
 gradT_excess_beta2 = s% gradT_excess_beta2
 gradT_excess_dlambda = s% gradT_excess_dlambda
 gradT_excess_dbeta = s% gradT_excess_dbeta
 
 D_mix_zero_region_bottom_q = s% D_mix_zero_region_bottom_q
 D_mix_zero_region_top_q = s% D_mix_zero_region_top_q
 dq_D_mix_zero_at_H_He_crossover = s% dq_D_mix_zero_at_H_He_crossover
 dq_D_mix_zero_at_H_C_crossover = s% dq_D_mix_zero_at_H_C_crossover

 use_superad_reduction = s% use_superad_reduction
 superad_reduction_gamma_limit = s% superad_reduction_gamma_limit
 superad_reduction_gamma_limit_scale = s% superad_reduction_gamma_limit_scale
 superad_reduction_gamma_inv_scale = s% superad_reduction_gamma_inv_scale
 superad_reduction_diff_grads_limit = s% superad_reduction_diff_grads_limit
 superad_reduction_limit = s% superad_reduction_limit
 
 max_logT_for_mlt = s% max_logT_for_mlt
 mlt_make_surface_no_mixing = s% mlt_make_surface_no_mixing
 do_normalize_dqs_as_part_of_set_qs = s% do_normalize_dqs_as_part_of_set_qs

 thermohaline_coeff = s% thermohaline_coeff
 thermohaline_option = s% thermohaline_option
 mixing_length_alpha = s% mixing_length_alpha
 remove_small_D_limit = s% remove_small_D_limit
 alt_scale_height_flag = s% alt_scale_height_flag
 Henyey_MLT_y_param = s% Henyey_MLT_y_param
 Henyey_MLT_nu_param = s% Henyey_MLT_nu_param
 make_gradr_sticky_in_solver_iters = s% make_gradr_sticky_in_solver_iters
 min_logT_for_make_gradr_sticky_in_solver_iters = s% min_logT_for_make_gradr_sticky_in_solver_iters
 no_MLT_below_shock = s% no_MLT_below_shock
 MLT_option = s% MLT_option
 mlt_use_rotation_correction = s% mlt_use_rotation_correction
 mlt_Pturb_factor = s% mlt_Pturb_factor

 burn_z_mix_region_logT = s% burn_z_mix_region_logT
 burn_he_mix_region_logT = s% burn_he_mix_region_logT
 burn_h_mix_region_logT = s% burn_h_mix_region_logT
 max_Y_for_burn_z_mix_region = s% max_Y_for_burn_z_mix_region
 max_X_for_burn_he_mix_region = s% max_X_for_burn_he_mix_region
 
 limit_overshoot_Hp_using_size_of_convection_zone = s% limit_overshoot_Hp_using_size_of_convection_zone

 predictive_mix = s%predictive_mix
 predictive_superad_thresh = s%predictive_superad_thresh
 predictive_avoid_reversal = s%predictive_avoid_reversal
 predictive_limit_ingestion = s%predictive_limit_ingestion
 predictive_ingestion_factor = s%predictive_ingestion_factor
 predictive_zone_type = s%predictive_zone_type
 predictive_zone_loc = s%predictive_zone_loc
 predictive_bdy_loc = s%predictive_bdy_loc
 predictive_bdy_q_min = s%predictive_bdy_q_min
 predictive_bdy_q_max = s%predictive_bdy_q_max

 do_conv_premix = s%do_conv_premix
 conv_premix_avoid_increase = s%conv_premix_avoid_increase
 conv_premix_time_factor = s%conv_premix_time_factor
 conv_premix_fix_pgas = s%conv_premix_fix_pgas
 conv_premix_dump_snapshots = s%conv_premix_dump_snapshots
 do_premix_heating = s%do_premix_heating

 overshoot_f = s%overshoot_f
 overshoot_f0 = s%overshoot_f0
 overshoot_D0 = s%overshoot_D0
 overshoot_Delta0 = s%overshoot_Delta0
 overshoot_mass_full_on = s%overshoot_mass_full_on
 overshoot_mass_full_off = s%overshoot_mass_full_off
 overshoot_scheme = s%overshoot_scheme
 overshoot_zone_type = s%overshoot_zone_type
 overshoot_zone_loc = s%overshoot_zone_loc
 overshoot_bdy_loc = s%overshoot_bdy_loc
 overshoot_D_min = s%overshoot_D_min
 overshoot_brunt_B_max = s%overshoot_brunt_B_max

 max_conv_vel_div_csound = s% max_conv_vel_div_csound
 max_v_for_convection = s% max_v_for_convection
 max_q_for_convection_with_hydro_on = s% max_q_for_convection_with_hydro_on
 max_v_div_cs_for_convection = s% max_v_div_cs_for_convection
 max_abs_du_div_cs_for_convection = s% max_abs_du_div_cs_for_convection

 calculate_Brunt_B = s% calculate_Brunt_B
 calculate_Brunt_N2 = s% calculate_Brunt_N2
 brunt_N2_coefficient = s% brunt_N2_coefficient
 threshold_for_smooth_brunt_B = s% threshold_for_smooth_brunt_B
 min_magnitude_brunt_B = s% min_magnitude_brunt_B

 min_overshoot_q = s% min_overshoot_q
 overshoot_alpha = s% overshoot_alpha
 
   RSP_max_num_periods = s% RSP_max_num_periods
   RSP_target_steps_per_cycle = s% RSP_target_steps_per_cycle
   RSP_min_max_R_for_periods = s% RSP_min_max_R_for_periods
   RSP_min_deltaR_for_periods = s% RSP_min_deltaR_for_periods
   RSP_default_PERIODLIN = s% RSP_default_PERIODLIN
   RSP_min_PERIOD_div_PERIODLIN = s% RSP_min_PERIOD_div_PERIODLIN
   RSP_GREKM_avg_abs_frac_new = s% RSP_GREKM_avg_abs_frac_new
   RSP_GREKM_avg_abs_limit = s% RSP_GREKM_avg_abs_limit
   RSP_theta = s% RSP_theta
   RSP_thetat = s% RSP_thetat
   RSP_thetau = s% RSP_thetau
   RSP_thetae = s% RSP_thetae
   RSP_thetaq = s% RSP_thetaq
   RSP_wtr = s% RSP_wtr
   RSP_wtc = s% RSP_wtc
   RSP_wtt = s% RSP_wtt
   RSP_gam = s% RSP_gam
   RSP_alfa = s% RSP_alfa
   RSP_alfap = s% RSP_alfap
   RSP_alfam = s% RSP_alfam
   RSP_alfat = s% RSP_alfat
   RSP_alfas = s% RSP_alfas
   RSP_alfac = s% RSP_alfac
   RSP_alfad = s% RSP_alfad
   RSP_gammar = s% RSP_gammar
   RSP_efl0 = s% RSP_efl0
   RSP_min_tau_for_turbulent_flux = s% RSP_min_tau_for_turbulent_flux
   RSP_cq = s% RSP_cq
   RSP_zsh = s% RSP_zsh
   RSP_Qvisc_quadratic = s% RSP_Qvisc_quadratic
   RSP_Qvisc_linear = s% RSP_Qvisc_linear
   RSP_Qvisc_linear_static = s% RSP_Qvisc_linear_static
   RSP_tol_max_corr = s% RSP_tol_max_corr
   RSP_tol_max_resid = s% RSP_tol_max_resid
   RSP_max_iters_per_try = s% RSP_max_iters_per_try
   RSP_max_retries_per_step = s% RSP_max_retries_per_step
   RSP_nz_div_IBOTOM = s% RSP_nz_div_IBOTOM
   RSP_kick_vsurf_km_per_sec = s% RSP_kick_vsurf_km_per_sec
   RSP_fraction_1st_overtone = s% RSP_fraction_1st_overtone
   RSP_fraction_2nd_overtone = s% RSP_fraction_2nd_overtone
   RSP_Avel = s% RSP_Avel
   RSP_Arnd = s% RSP_Arnd
   RSP_mode_for_setting_PERIODLIN = s% RSP_mode_for_setting_PERIODLIN
   RSP_initial_dt_factor = s% RSP_initial_dt_factor
   RSP_v_div_cs_threshold_for_dt_limit = s% RSP_v_div_cs_threshold_for_dt_limit
   RSP_max_dt_times_min_dr_div_cs = s% RSP_max_dt_times_min_dr_div_cs
   RSP_max_dt_times_min_rad_diff_time = s% RSP_max_dt_times_min_rad_diff_time
   RSP_max_dt = s% RSP_max_dt
   RSP_testing = s% RSP_testing
   RSP_report_limit_dt = s% RSP_report_limit_dt
   RSP_use_Prad_for_Psurf = s% RSP_use_Prad_for_Psurf
   RSP_report_undercorrections = s% RSP_report_undercorrections
   RSP_use_atm_grey_with_kap_for_Psurf = s% RSP_use_atm_grey_with_kap_for_Psurf
   use_other_RSP_linear_analysis = s% use_other_RSP_linear_analysis
   use_other_RSP_build_model = s% use_other_RSP_build_model
   RSP_kap_density_factor = s% RSP_kap_density_factor
   RSP_fixed_Psurf = s% RSP_fixed_Psurf
   RSP_hydro_only = s% RSP_hydro_only
   RSP_tau_surf_for_atm_grey_with_kap = s% RSP_tau_surf_for_atm_grey_with_kap
   RSP_Psurf = s% RSP_Psurf
   set_RSP_Psurf_to_multiple_of_initial_P1 = s% set_RSP_Psurf_to_multiple_of_initial_P1
   RSP_surface_tau = s% RSP_surface_tau
   RSP_write_map = s% RSP_write_map
   RSP_trace_RSP_build_model = s% RSP_trace_RSP_build_model
   RSP_map_filename = s% RSP_map_filename
   RSP_map_columns_filename = s% RSP_map_columns_filename
   RSP_map_history_filename = s% RSP_map_history_filename
   RSP_map_first_period = s% RSP_map_first_period
   RSP_map_last_period = s% RSP_map_last_period
   RSP_map_zone_interval = s% RSP_map_zone_interval
   RSP_nmodes = s% RSP_nmodes
   RSP_work_period = s% RSP_work_period
   RSP_work_filename = s% RSP_work_filename
   RSP_nz_outer = s% RSP_nz_outer
   RSP_max_outer_dm_tries = s% RSP_max_outer_dm_tries
   RSP_max_inner_scale_tries = s% RSP_max_inner_scale_tries
   RSP_relax_max_tries = s% RSP_relax_max_tries
   RSP_T_anchor_tolerance = s% RSP_T_anchor_tolerance
   RSP_T_inner_tolerance = s% RSP_T_inner_tolerance
   RSP_relax_dm_tolerance = s% RSP_relax_dm_tolerance
   RSP_dq_1_factor = s% RSP_dq_1_factor
   use_RSP_new_start_scheme = s% use_RSP_new_start_scheme
   RSP_do_check_omega = s% RSP_do_check_omega
   RSP_report_check_omega_changes = s% RSP_report_check_omega_changes
   RSP_nz = s% RSP_nz
   RSP_T_anchor = s% RSP_T_anchor
   RSP_T_inner = s% RSP_T_inner
   RSP_relax_initial_model = s% RSP_relax_initial_model
   RSP_relax_alfap_before_alfat = s% RSP_relax_alfap_before_alfat
   RSP_relax_adjust_inner_mass_distribution = s% RSP_relax_adjust_inner_mass_distribution
   RSP_Teff = s% RSP_Teff
   RSP_mass = s% RSP_mass
   RSP_L = s% RSP_L
   RSP_X = s% RSP_X
   RSP_Z = s% RSP_Z

 RTI_smooth_mass = s% RTI_smooth_mass
 RTI_smooth_iterations = s% RTI_smooth_iterations
 RTI_smooth_fraction = s% RTI_smooth_fraction

 alpha_RTI_diffusion_factor = s% alpha_RTI_diffusion_factor
 dudt_RTI_diffusion_factor = s% dudt_RTI_diffusion_factor
 dedt_RTI_diffusion_factor = s% dedt_RTI_diffusion_factor
 dlnddt_RTI_diffusion_factor = s% dlnddt_RTI_diffusion_factor
 composition_RTI_diffusion_factor = s% composition_RTI_diffusion_factor
 max_M_RTI_factors_full_on = s% max_M_RTI_factors_full_on
 min_M_RTI_factors_full_off = s% min_M_RTI_factors_full_off

 alpha_RTI_src_min_v_div_cs = s% alpha_RTI_src_min_v_div_cs
 alpha_RTI_src_max_q = s% alpha_RTI_src_max_q
 alpha_RTI_src_min_q = s% alpha_RTI_src_min_q

 T_mix_limit = s% T_mix_limit
 mlt_gradT_fraction = s% mlt_gradT_fraction

 ! atmosphere -- surface boundary conditions
 atm_option = s% atm_option
 atm_off_table_option = s% atm_off_table_option
 Pextra_factor = s% Pextra_factor
 atm_fixed_Teff = s% atm_fixed_Teff
 atm_fixed_Psurf = s% atm_fixed_Psurf
 atm_fixed_Tsurf = s% atm_fixed_Tsurf

 atm_T_tau_relation = s% atm_T_tau_relation
 atm_T_tau_opacity = s% atm_T_tau_opacity
 atm_T_tau_errtol = s% atm_T_tau_errtol
 atm_T_tau_max_iters = s% atm_T_tau_max_iters
 atm_T_tau_max_steps = s% atm_T_tau_max_steps

 atm_table = s% atm_table

 atm_irradiated_opacity = s% atm_irradiated_opacity
 atm_irradiated_errtol = s% atm_irradiated_errtol
 atm_irradiated_T_eq = s% atm_irradiated_T_eq
 atm_irradiated_kap_v = s% atm_irradiated_kap_v
 atm_irradiated_kap_v_div_kap_th = s% atm_irradiated_kap_v_div_kap_th
 atm_irradiated_P_surf = s% atm_irradiated_P_surf
 atm_irradiated_max_iters = s% atm_irradiated_max_iters

 use_compression_outer_BC = s% use_compression_outer_BC
 use_momentum_outer_BC = s% use_momentum_outer_BC
 Tsurf_factor = s% Tsurf_factor
 use_zero_Pgas_outer_BC = s% use_zero_Pgas_outer_BC
 fixed_vsurf = s% fixed_vsurf
 use_fixed_vsurf_outer_BC = s% use_fixed_vsurf_outer_BC
 fixed_Psurf = s% fixed_Psurf
 use_fixed_Psurf_outer_BC = s% use_fixed_Psurf_outer_BC

 atm_build_tau_outer = s% atm_build_tau_outer
 atm_build_dlogtau = s% atm_build_dlogtau
 atm_build_errtol = s% atm_build_errtol
 
 use_T_tau_gradr_factor = s% use_T_tau_gradr_factor

 ! extra heat near surface to model irradiation
 irradiation_flux = s% irradiation_flux
 column_depth_for_irradiation = s% column_depth_for_irradiation

 ! extra heat
 inject_uniform_extra_heat = s% inject_uniform_extra_heat
 min_q_for_uniform_extra_heat = s% min_q_for_uniform_extra_heat
 max_q_for_uniform_extra_heat = s% max_q_for_uniform_extra_heat
 inject_extra_ergs_sec = s% inject_extra_ergs_sec
 base_of_inject_extra_ergs_sec = s% base_of_inject_extra_ergs_sec
 total_mass_for_inject_extra_ergs_sec = s% total_mass_for_inject_extra_ergs_sec
 start_time_for_inject_extra_ergs_sec = s% start_time_for_inject_extra_ergs_sec
 duration_for_inject_extra_ergs_sec = s% duration_for_inject_extra_ergs_sec
 inject_until_reach_model_with_total_energy = s% inject_until_reach_model_with_total_energy

 ! mass gain or loss
 mass_change = s% mass_change
 mass_change_full_off_dt = s% mass_change_full_off_dt
 mass_change_full_on_dt = s% mass_change_full_on_dt
 trace_dt_control_mass_change = s% trace_dt_control_mass_change
 no_wind_if_no_rotation = s% no_wind_if_no_rotation

 min_wind = s% min_wind
 max_wind = s% max_wind
 use_accreted_material_j = s% use_accreted_material_j
 accreted_material_j = s% accreted_material_j
 D_omega_mixing_rate = s% D_omega_mixing_rate
 D_omega_mixing_across_convection_boundary = s% D_omega_mixing_across_convection_boundary
 max_q_for_D_omega_zero_in_convection_region = s% max_q_for_D_omega_zero_in_convection_region
 nu_omega_mixing_rate = s% nu_omega_mixing_rate
 nu_omega_mixing_across_convection_boundary = s% nu_omega_mixing_across_convection_boundary
 max_q_for_nu_omega_zero_in_convection_region = s% max_q_for_nu_omega_zero_in_convection_region
 
 mdot_omega_power = s% mdot_omega_power
 max_rotational_mdot_boost = s% max_rotational_mdot_boost
 max_mdot_jump_for_rotation = s% max_mdot_jump_for_rotation
 lim_trace_rotational_mdot_boost = s% lim_trace_rotational_mdot_boost
 rotational_mdot_boost_fac = s% rotational_mdot_boost_fac
 rotational_mdot_kh_fac = s% rotational_mdot_kh_fac
 surf_avg_tau = s% surf_avg_tau
 surf_avg_tau_min = s% surf_avg_tau_min

 super_eddington_scaling_factor = s% super_eddington_scaling_factor
 super_eddington_wind_Ledd_factor = s% super_eddington_wind_Ledd_factor
 wind_boost_full_off_L_div_Ledd = s% wind_boost_full_off_L_div_Ledd
 wind_boost_full_on_L_div_Ledd = s% wind_boost_full_on_L_div_Ledd
 super_eddington_wind_max_boost = s% super_eddington_wind_max_boost
 trace_super_eddington_wind_boost = s% trace_super_eddington_wind_boost
 
 max_tries_for_implicit_wind = s% max_tries_for_implicit_wind
 iwind_tolerance = s% iwind_tolerance
 iwind_lambda = s% iwind_lambda

 rlo_scaling_factor = s% rlo_scaling_factor
 rlo_wind_min_L = s% rlo_wind_min_L
 rlo_wind_max_Teff = s% rlo_wind_max_Teff
 rlo_wind_roche_lobe_radius = s% rlo_wind_roche_lobe_radius
 roche_lobe_xfer_full_on = s% roche_lobe_xfer_full_on
 roche_lobe_xfer_full_off = s% roche_lobe_xfer_full_off
 rlo_wind_base_mdot = s% rlo_wind_base_mdot
 rlo_wind_scale_height = s% rlo_wind_scale_height

 cool_wind_RGB_scheme = s% cool_wind_RGB_scheme
 cool_wind_AGB_scheme = s% cool_wind_AGB_scheme
 RGB_to_AGB_wind_switch = s% RGB_to_AGB_wind_switch
 Reimers_scaling_factor = s% Reimers_scaling_factor
 Blocker_scaling_factor = s% Blocker_scaling_factor
 de_Jager_scaling_factor = s% de_Jager_scaling_factor
 van_Loon_scaling_factor = s% van_Loon_scaling_factor
 Nieuwenhuijzen_scaling_factor = s% Nieuwenhuijzen_scaling_factor
 Vink_scaling_factor = s% Vink_scaling_factor
 Dutch_scaling_factor = s% Dutch_scaling_factor
 Dutch_wind_lowT_scheme = s% Dutch_wind_lowT_scheme

 wind_H_envelope_limit = s% wind_H_envelope_limit
 wind_H_He_envelope_limit = s% wind_H_He_envelope_limit
 wind_He_layer_limit = s% wind_He_layer_limit

 max_logT_for_k_below_const_q = s% max_logT_for_k_below_const_q
 max_q_for_k_below_const_q = s% max_q_for_k_below_const_q
 min_q_for_k_below_const_q = s% min_q_for_k_below_const_q
 max_logT_for_k_const_mass = s% max_logT_for_k_const_mass
 min_q_for_k_const_mass = s% min_q_for_k_const_mass
 max_q_for_k_const_mass = s% max_q_for_k_const_mass

 ! composition of added mass
 accrete_same_as_surface = s% accrete_same_as_surface

 accrete_given_mass_fractions = s% accrete_given_mass_fractions
 num_accretion_species = s% num_accretion_species
 accretion_species_id = s% accretion_species_id
 accretion_species_xa = s% accretion_species_xa

 accretion_h1 = s% accretion_h1
 accretion_h2 = s% accretion_h2
 accretion_he3 = s% accretion_he3
 accretion_he4 = s% accretion_he4
 accretion_zfracs = s% accretion_zfracs
 accretion_dump_missing_metals_into_heaviest = s% accretion_dump_missing_metals_into_heaviest

 ! special list of z fractions
 z_fraction_li = s% z_fraction_li
 z_fraction_be = s% z_fraction_be
 z_fraction_b = s% z_fraction_b
 z_fraction_c = s% z_fraction_c
 z_fraction_n = s% z_fraction_n
 z_fraction_o = s% z_fraction_o
 z_fraction_f = s% z_fraction_f
 z_fraction_ne = s% z_fraction_ne
 z_fraction_na = s% z_fraction_na
 z_fraction_mg = s% z_fraction_mg
 z_fraction_al = s% z_fraction_al
 z_fraction_si = s% z_fraction_si
 z_fraction_p = s% z_fraction_p
 z_fraction_s = s% z_fraction_s
 z_fraction_cl = s% z_fraction_cl
 z_fraction_ar = s% z_fraction_ar
 z_fraction_k = s% z_fraction_k
 z_fraction_ca = s% z_fraction_ca
 z_fraction_sc = s% z_fraction_sc
 z_fraction_ti = s% z_fraction_ti
 z_fraction_v = s% z_fraction_v
 z_fraction_cr = s% z_fraction_cr
 z_fraction_mn = s% z_fraction_mn
 z_fraction_fe = s% z_fraction_fe
 z_fraction_co = s% z_fraction_co
 z_fraction_ni = s% z_fraction_ni
 z_fraction_cu = s% z_fraction_cu
 z_fraction_zn = s% z_fraction_zn

 lgT_lo_for_set_new_abundances = s% lgT_lo_for_set_new_abundances
 lgT_hi_for_set_new_abundances = s% lgT_hi_for_set_new_abundances

 ! automatic stops for mass loss/gain
 max_star_mass_for_gain = s% max_star_mass_for_gain
 min_star_mass_for_loss = s% min_star_mass_for_loss
 max_T_center_for_any_mass_loss = s% max_T_center_for_any_mass_loss
 max_T_center_for_full_mass_loss = s% max_T_center_for_full_mass_loss

 ! relaxation parameters
 extra_power_source = s% extra_power_source
 relax_dlnZ = s% relax_dlnZ
 relax_dY = s% relax_dY

 ! mesh adjustment
 show_mesh_changes = s% show_mesh_changes
 okay_to_remesh = s% okay_to_remesh
 restore_mesh_on_retry = s% restore_mesh_on_retry
 num_steps_to_hold_mesh_after_retry = s% num_steps_to_hold_mesh_after_retry
 trace_mesh_adjust_error_in_conservation = s% trace_mesh_adjust_error_in_conservation
 max_rel_delta_IE_for_mesh_total_energy_balance = s% max_rel_delta_IE_for_mesh_total_energy_balance
 max_allowed_nz = s% max_allowed_nz
 mesh_max_allowed_ratio = s% mesh_max_allowed_ratio
 remesh_max_allowed_logT = s% remesh_max_allowed_logT
 max_delta_x_for_merge = s% max_delta_x_for_merge

 mesh_ok_to_merge = s% mesh_ok_to_merge
 mesh_max_k_old_for_split = s% mesh_max_k_old_for_split
 mesh_min_k_old_for_split = s% mesh_min_k_old_for_split
 mesh_adjust_get_T_from_E = s% mesh_adjust_get_T_from_E

 max_dq = s% max_dq
 min_dq = s% min_dq
 min_dq_for_split = s% min_dq_for_split
 min_dq_for_xa = s% min_dq_for_xa
 min_dq_for_xa_convective = s% min_dq_for_xa_convective
 min_dq_for_logT = s% min_dq_for_logT

 mesh_min_dlnR = s% mesh_min_dlnR
 merge_if_dlnR_too_small = s% merge_if_dlnR_too_small

 mesh_min_dr_div_dRstar = s% mesh_min_dr_div_dRstar
 merge_if_dr_div_dRstar_too_small = s% merge_if_dr_div_dRstar_too_small

 mesh_min_dr_div_cs = s% mesh_min_dr_div_cs
 merge_if_dr_div_cs_too_small = s% merge_if_dr_div_cs_too_small

 max_center_cell_dq = s% max_center_cell_dq
 max_surface_cell_dq = s% max_surface_cell_dq
 max_num_subcells = s% max_num_subcells
 max_num_merge_cells = s% max_num_merge_cells

 mesh_delta_coeff = s% mesh_delta_coeff
 mesh_delta_coeff_for_highT = s% mesh_delta_coeff_for_highT
 logT_max_for_standard_mesh_delta_coeff = s% logT_max_for_standard_mesh_delta_coeff
 logT_min_for_highT_mesh_delta_coeff = s% logT_min_for_highT_mesh_delta_coeff
 mesh_Pgas_div_P_exponent = s% mesh_Pgas_div_P_exponent

 remesh_dt_limit = s% remesh_dt_limit

 E_function_weight = s% E_function_weight
 E_function_param = s% E_function_param
 P_function_weight = s% P_function_weight

 mesh_logX_species = s% mesh_logX_species
 mesh_logX_min_for_extra = s% mesh_logX_min_for_extra
 mesh_dlogX_dlogP_extra = s% mesh_dlogX_dlogP_extra
 mesh_dlogX_dlogP_full_on = s% mesh_dlogX_dlogP_full_on
 mesh_dlogX_dlogP_full_off = s% mesh_dlogX_dlogP_full_off

 mesh_dlog_eps_min_for_extra = s% mesh_dlog_eps_min_for_extra
 mesh_dlog_eps_dlogP_full_on = s% mesh_dlog_eps_dlogP_full_on
 mesh_dlog_eps_dlogP_full_off = s% mesh_dlog_eps_dlogP_full_off

 mesh_dlog_pp_dlogP_extra = s% mesh_dlog_pp_dlogP_extra
 mesh_dlog_cno_dlogP_extra = s% mesh_dlog_cno_dlogP_extra
 mesh_dlog_3alf_dlogP_extra = s% mesh_dlog_3alf_dlogP_extra

 mesh_dlog_burn_c_dlogP_extra = s% mesh_dlog_burn_c_dlogP_extra
 mesh_dlog_burn_n_dlogP_extra = s% mesh_dlog_burn_n_dlogP_extra
 mesh_dlog_burn_o_dlogP_extra = s% mesh_dlog_burn_o_dlogP_extra
 mesh_dlog_burn_ne_dlogP_extra = s% mesh_dlog_burn_ne_dlogP_extra
 mesh_dlog_burn_na_dlogP_extra = s% mesh_dlog_burn_na_dlogP_extra
 mesh_dlog_burn_mg_dlogP_extra = s% mesh_dlog_burn_mg_dlogP_extra
 mesh_dlog_burn_si_dlogP_extra = s% mesh_dlog_burn_si_dlogP_extra
 mesh_dlog_burn_s_dlogP_extra = s% mesh_dlog_burn_s_dlogP_extra
 mesh_dlog_burn_ar_dlogP_extra = s% mesh_dlog_burn_ar_dlogP_extra
 mesh_dlog_burn_ca_dlogP_extra = s% mesh_dlog_burn_ca_dlogP_extra
 mesh_dlog_burn_ti_dlogP_extra = s% mesh_dlog_burn_ti_dlogP_extra
 mesh_dlog_burn_cr_dlogP_extra = s% mesh_dlog_burn_cr_dlogP_extra
 mesh_dlog_burn_fe_dlogP_extra = s% mesh_dlog_burn_fe_dlogP_extra

 mesh_dlog_cc_dlogP_extra = s% mesh_dlog_cc_dlogP_extra
 mesh_dlog_co_dlogP_extra = s% mesh_dlog_co_dlogP_extra
 mesh_dlog_oo_dlogP_extra = s% mesh_dlog_oo_dlogP_extra

 mesh_dlog_pnhe4_dlogP_extra = s% mesh_dlog_pnhe4_dlogP_extra
 mesh_dlog_photo_dlogP_extra = s% mesh_dlog_photo_dlogP_extra
 mesh_dlog_other_dlogP_extra = s% mesh_dlog_other_dlogP_extra
 
 mesh_delta_coeff_factor_smooth_iters = s% mesh_delta_coeff_factor_smooth_iters

 T_function1_weight = s% T_function1_weight
 T_function2_weight = s% T_function2_weight
 T_function2_param = s% T_function2_param

 R_function_weight = s% R_function_weight
 R_function_param = s% R_function_param

 R_function2_weight = s% R_function2_weight
 R_function2_param1 = s% R_function2_param1
 R_function2_param2 = s% R_function2_param2

 R_function3_weight = s% R_function3_weight

 M_function_weight = s% M_function_weight
 M_function_param = s% M_function_param

 gradT_function_weight = s% gradT_function_weight
 log_tau_function_weight = s% log_tau_function_weight
 log_kap_function_weight = s% log_kap_function_weight
 omega_function_weight = s% omega_function_weight

 gam_function_weight = s% gam_function_weight
 gam_function_param1 = s% gam_function_param1
 gam_function_param2 = s% gam_function_param2

 xa_function_species = s% xa_function_species
 xa_function_weight = s% xa_function_weight
 xa_function_param = s% xa_function_param
 xa_mesh_delta_coeff = s% xa_mesh_delta_coeff
 
 use_split_merge_amr = s% use_split_merge_amr
 split_merge_amr_nz_baseline = s% split_merge_amr_nz_baseline
 split_merge_amr_nz_r_core = s% split_merge_amr_nz_r_core
 split_merge_amr_nz_r_core_fraction = s% split_merge_amr_nz_r_core_fraction
 split_merge_amr_mesh_delta_coeff = s% split_merge_amr_mesh_delta_coeff
 split_merge_amr_log_zoning = s% split_merge_amr_log_zoning
 split_merge_amr_hybrid_zoning = s% split_merge_amr_hybrid_zoning
 split_merge_amr_flipped_hybrid_zoning = s% split_merge_amr_flipped_hybrid_zoning
 split_merge_amr_logtau_zoning = s% split_merge_amr_logtau_zoning
 split_merge_amr_okay_to_split_nz = s% split_merge_amr_okay_to_split_nz
 split_merge_amr_okay_to_split_1 = s% split_merge_amr_okay_to_split_1
 merge_amr_inhibit_at_jumps = s% merge_amr_inhibit_at_jumps
 split_merge_amr_MaxLong = s% split_merge_amr_MaxLong
 split_merge_amr_MaxShort = s% split_merge_amr_MaxShort
 merge_amr_max_abs_du_div_cs = s% merge_amr_max_abs_du_div_cs
 merge_amr_ignore_surface_cells = s% merge_amr_ignore_surface_cells
 merge_amr_du_div_cs_limit_only_for_compression = s% merge_amr_du_div_cs_limit_only_for_compression
 split_merge_amr_avoid_repeated_remesh = s% split_merge_amr_avoid_repeated_remesh
 merge_amr_k_for_ignore_surface_cells = s% merge_amr_k_for_ignore_surface_cells
 split_merge_amr_dq_min = s% split_merge_amr_dq_min
 split_merge_amr_dq_max = s% split_merge_amr_dq_max
 split_merge_amr_r_core_cm = s% split_merge_amr_r_core_cm
 split_merge_amr_max_iters = s% split_merge_amr_max_iters
 trace_split_merge_amr = s% trace_split_merge_amr
 equal_split_density_amr = s% equal_split_density_amr

 ! nuclear reaction parameters
 screening_mode = s% screening_mode
 default_net_name = s% default_net_name

 net_logTcut_lo = s% net_logTcut_lo
 net_logTcut_lim = s% net_logTcut_lim

 eps_nuc_factor = s% eps_nuc_factor
 op_split_burn_eps_nuc_infall_limit = s% op_split_burn_eps_nuc_infall_limit
 eps_WD_sedimentation_factor = s% eps_WD_sedimentation_factor
 max_abs_eps_nuc = s% max_abs_eps_nuc
 dxdt_nuc_factor = s% dxdt_nuc_factor
 max_abar_for_burning = s% max_abar_for_burning
 fe56ec_fake_factor = s% fe56ec_fake_factor
 min_T_for_fe56ec_fake_factor = s% min_T_for_fe56ec_fake_factor
 weak_rate_factor = s% weak_rate_factor

 mix_factor = s% mix_factor

 sig_term_limit = s% sig_term_limit

 sig_min_factor_for_high_Tcenter = s% sig_min_factor_for_high_Tcenter
 Tcenter_min_for_sig_min_factor_full_on = s% Tcenter_min_for_sig_min_factor_full_on
 Tcenter_max_for_sig_min_factor_full_off = s% Tcenter_max_for_sig_min_factor_full_off
 max_delta_m_to_bdy_for_sig_min_factor = s% max_delta_m_to_bdy_for_sig_min_factor
 delta_m_lower_for_sig_min_factor = s% delta_m_lower_for_sig_min_factor
 delta_m_upper_for_sig_min_factor = s% delta_m_upper_for_sig_min_factor

 am_sig_term_limit = s% am_sig_term_limit
 am_D_mix_factor = s% am_D_mix_factor
 am_gradmu_factor = s% am_gradmu_factor
 am_nu_factor = s% am_nu_factor

 D_visc_factor = s% D_visc_factor
 D_DSI_factor = s% D_DSI_factor
 D_SH_factor = s% D_SH_factor
 D_SSI_factor = s% D_SSI_factor
 D_ES_factor = s% D_ES_factor
 D_GSF_factor = s% D_GSF_factor
 D_ST_factor = s% D_ST_factor

 am_nu_non_rotation_factor = s% am_nu_non_rotation_factor
 skip_rotation_in_convection_zones = s% skip_rotation_in_convection_zones
 am_nu_DSI_factor = s% am_nu_DSI_factor
 am_nu_SH_factor = s% am_nu_SH_factor
 am_nu_SSI_factor = s% am_nu_SSI_factor
 am_nu_ES_factor = s% am_nu_ES_factor
 am_nu_GSF_factor = s% am_nu_GSF_factor
 am_nu_ST_factor = s% am_nu_ST_factor
 am_nu_visc_factor = s% am_nu_visc_factor

 am_nu_omega_rot_factor = s% am_nu_omega_rot_factor
 am_nu_omega_non_rot_factor = s% am_nu_omega_non_rot_factor
 am_nu_j_rot_factor = s% am_nu_j_rot_factor
 am_nu_j_non_rot_factor = s% am_nu_j_non_rot_factor

 smooth_nu_ST = s% smooth_nu_ST
 smooth_D_ST = s% smooth_D_ST
 smooth_D_DSI = s% smooth_D_DSI
 smooth_D_SSI = s% smooth_D_SSI
 smooth_D_SH = s% smooth_D_SH
 smooth_D_GSF = s% smooth_D_GSF
 smooth_D_ES = s% smooth_D_ES
 smooth_D_omega = s% smooth_D_omega
 smooth_am_nu_rot = s% smooth_am_nu_rot

 simple_i_rot_flag = s% simple_i_rot_flag
 do_adjust_J_lost = s% do_adjust_J_lost
 premix_omega = s% premix_omega
 angular_momentum_error_warn = s% angular_momentum_error_warn
 angular_momentum_error_retry = s% angular_momentum_error_retry
 recalc_mixing_info_each_substep = s% recalc_mixing_info_each_substep
 adjust_J_fraction = s% adjust_J_fraction
 min_q_for_adjust_J_lost = s% min_q_for_adjust_J_lost
 min_J_div_delta_J = s% min_J_div_delta_J
 max_mdot_redo_cnt = s% max_mdot_redo_cnt
 mdot_revise_factor = s% mdot_revise_factor
 implicit_mdot_boost = s% implicit_mdot_boost
 min_years_dt_for_redo_mdot = s% min_years_dt_for_redo_mdot
 surf_omega_div_omega_crit_limit = s% surf_omega_div_omega_crit_limit
 surf_omega_div_omega_crit_tol = S% surf_omega_div_omega_crit_tol
 fitted_fp_ft_i_rot = s% fitted_fp_ft_i_rot 
 w_div_wcrit_max = s% w_div_wcrit_max
 w_div_wcrit_max2 = s% w_div_wcrit_max2
 fp_min = s% fp_min
 ft_min = s% ft_min
 fp_error_limit = s% fp_error_limit
 ft_error_limit = s% ft_error_limit

 D_mix_rotation_max_logT_full_on = s% D_mix_rotation_max_logT_full_on
 D_mix_rotation_min_logT_full_off = s% D_mix_rotation_min_logT_full_off

 set_uniform_am_nu_non_rot = s% set_uniform_am_nu_non_rot
 uniform_am_nu_non_rot = s% uniform_am_nu_non_rot

 set_min_am_nu_non_rot = s% set_min_am_nu_non_rot
 min_am_nu_non_rot = s% min_am_nu_non_rot
 min_center_Ye_for_min_am_nu_non_rot = s% min_center_Ye_for_min_am_nu_non_rot

 set_min_D_mix = s% set_min_D_mix
 mass_lower_limit_for_min_D_mix = s% mass_lower_limit_for_min_D_mix
 mass_upper_limit_for_min_D_mix = s% mass_upper_limit_for_min_D_mix
 min_D_mix = s% min_D_mix
 set_min_D_mix_below_Tmax = s% set_min_D_mix_below_Tmax
 min_D_mix_below_Tmax = s% min_D_mix_below_Tmax
 set_min_D_mix_in_H_He = s% set_min_D_mix_in_H_He
 min_D_mix_in_H_He = s% min_D_mix_in_H_He
 min_center_Ye_for_min_D_mix = s% min_center_Ye_for_min_D_mix
 reaction_neuQs_factor = s% reaction_neuQs_factor
 nonlocal_NiCo_kap_gamma = s% nonlocal_NiCo_kap_gamma
 nonlocal_NiCo_decay_heat = s% nonlocal_NiCo_decay_heat
 dtau_gamma_NiCo_decay_heat = s% dtau_gamma_NiCo_decay_heat
 max_logT_for_net = s% max_logT_for_net
 smooth_outer_xa_big = s% smooth_outer_xa_big
 smooth_outer_xa_small = s% smooth_outer_xa_small

 ! element diffusion parameters
 diffusion_use_iben_macdonald = s% diffusion_use_iben_macdonald
 diffusion_use_paquette = s% diffusion_use_paquette
 diffusion_use_cgs_solver = s% diffusion_use_cgs_solver
 diffusion_use_full_net = s% diffusion_use_full_net
 do_WD_sedimentation_heating = s% do_WD_sedimentation_heating
 min_xa_for_WD_sedimentation_heating = s% min_xa_for_WD_sedimentation_heating
 do_diffusion_heating = s% do_diffusion_heating
 do_element_diffusion = s% do_element_diffusion
 cgs_thermal_diffusion_eta_full_on = s% cgs_thermal_diffusion_eta_full_on
 cgs_thermal_diffusion_eta_full_off = s% cgs_thermal_diffusion_eta_full_off
 diffusion_min_dq_at_surface = s% diffusion_min_dq_at_surface
 diffusion_min_T_at_surface = s% diffusion_min_T_at_surface
 diffusion_min_dq_ratio_at_surface = s% diffusion_min_dq_ratio_at_surface
 diffusion_dt_limit = s% diffusion_dt_limit

 diffusion_min_X_hard_limit = s% diffusion_min_X_hard_limit
 diffusion_X_total_atol = s% diffusion_X_total_atol
 diffusion_X_total_rtol = s% diffusion_X_total_rtol
 diffusion_upwind_abs_v_limit = s% diffusion_upwind_abs_v_limit
 diffusion_dt_div_timescale = s% diffusion_dt_div_timescale
 diffusion_min_num_substeps = s% diffusion_min_num_substeps
 diffusion_max_iters_per_substep = s% diffusion_max_iters_per_substep
 diffusion_max_retries_per_substep = s% diffusion_max_retries_per_substep
 diffusion_v_max = s% diffusion_v_max
 diffusion_gamma_full_off = s% diffusion_gamma_full_off
 diffusion_gamma_full_on = s% diffusion_gamma_full_on
 diffusion_T_full_off = s% diffusion_T_full_off
 D_mix_ignore_diffusion = s% D_mix_ignore_diffusion
 diffusion_T_full_on = s% diffusion_T_full_on
 diffusion_calculates_ionization = s% diffusion_calculates_ionization
 diffusion_nsmooth_typical_charge = s% diffusion_nsmooth_typical_charge
 diffusion_tol_correction_max = s% diffusion_tol_correction_max
 diffusion_tol_correction_norm = s% diffusion_tol_correction_norm

 diffusion_AD_dm_full_on = s% diffusion_AD_dm_full_on
 diffusion_AD_dm_full_off = s% diffusion_AD_dm_full_off
 diffusion_AD_boost_factor = s% diffusion_AD_boost_factor

 diffusion_SIG_factor = s% diffusion_SIG_factor
 diffusion_GT_factor = s% diffusion_GT_factor

 diffusion_Vlimit_dm_full_on = s% diffusion_Vlimit_dm_full_on
 diffusion_Vlimit_dm_full_off = s% diffusion_Vlimit_dm_full_off
 diffusion_Vlimit = s% diffusion_Vlimit

 diffusion_max_T_for_radaccel = s% diffusion_max_T_for_radaccel
 diffusion_min_T_for_radaccel = s% diffusion_min_T_for_radaccel
 diffusion_max_Z_for_radaccel = s% diffusion_max_Z_for_radaccel
 diffusion_min_Z_for_radaccel = s% diffusion_min_Z_for_radaccel
 diffusion_screening_for_radaccel = s% diffusion_screening_for_radaccel
 op_mono_data_path = s% op_mono_data_path
 op_mono_data_cache_filename = s% op_mono_data_cache_filename

 show_diffusion_info = s% show_diffusion_info
 show_diffusion_substep_info = s% show_diffusion_substep_info
 show_diffusion_timing = s% show_diffusion_timing

 diffusion_num_classes = s% diffusion_num_classes
 diffusion_class_representative = s% diffusion_class_representative
 diffusion_class_A_max = s% diffusion_class_A_max
 diffusion_class_typical_charge = s% diffusion_class_typical_charge
 diffusion_class_factor = s% diffusion_class_factor

 diffusion_use_isolve = s% diffusion_use_isolve
 diffusion_rtol_for_isolve = s% diffusion_rtol_for_isolve
 diffusion_atol_for_isolve = s% diffusion_atol_for_isolve
 diffusion_maxsteps_for_isolve = s% diffusion_maxsteps_for_isolve
 diffusion_isolve_solver = s% diffusion_isolve_solver

 ! eos controls
 use_fixed_XZ_for_eos = s% use_fixed_XZ_for_eos
 fixed_X_for_eos = s% fixed_X_for_eos
 fixed_Z_for_eos = s% fixed_Z_for_eos
 report_eos_settings_at_start_of_run = s% report_eos_settings_at_start_of_run
 use_d_eos_dxa = s% use_d_eos_dxa
 
 ! opacity controls
 use_simple_es_for_kap = s% use_simple_es_for_kap
 use_starting_composition_for_kap = s% use_starting_composition_for_kap
 min_kap_for_dPrad_dm_eqn = s% min_kap_for_dPrad_dm_eqn

 low_logT_op_mono_full_off = s% low_logT_op_mono_full_off
 low_logT_op_mono_full_on = s% low_logT_op_mono_full_on
 high_logT_op_mono_full_off = s% high_logT_op_mono_full_off
 high_logT_op_mono_full_on = s% high_logT_op_mono_full_on
 op_mono_min_X_to_include = s% op_mono_min_X_to_include
 use_op_mono_alt_get_kap = s% use_op_mono_alt_get_kap

 include_L_in_correction_limits = s% include_L_in_correction_limits
 include_v_in_correction_limits = s% include_v_in_correction_limits
 include_u_in_correction_limits = s% include_u_in_correction_limits
 include_w_in_correction_limits = s% include_w_in_correction_limits

 ! asteroseismology controls

 get_delta_nu_from_scaled_solar = s% get_delta_nu_from_scaled_solar
 nu_max_sun = s% nu_max_sun
 delta_nu_sun = s% delta_nu_sun
 Teff_sun = s% Teff_sun
 delta_Pg_mode_freq = s% delta_Pg_mode_freq

 ! hydro parameters
 opacity_max = s% opacity_max
 opacity_factor = s% opacity_factor
 min_logT_for_opacity_factor_off = s% min_logT_for_opacity_factor_off
 min_logT_for_opacity_factor_on = s% min_logT_for_opacity_factor_on
 max_logT_for_opacity_factor_on = s% max_logT_for_opacity_factor_on
 max_logT_for_opacity_factor_off = s% max_logT_for_opacity_factor_off

 non_nuc_neu_factor = s% non_nuc_neu_factor
 use_dedt_form_of_energy_eqn = s% use_dedt_form_of_energy_eqn
 always_use_dedt_form_of_energy_eqn = s% always_use_dedt_form_of_energy_eqn
 use_time_centered_eps_grav = s% use_time_centered_eps_grav
 no_dedt_form_during_relax = s% no_dedt_form_during_relax
 always_use_eps_grav_form_of_energy_eqn = s% always_use_eps_grav_form_of_energy_eqn
 dedt_eqn_r_scale = s% dedt_eqn_r_scale
 use_mass_corrections = s% use_mass_corrections
 use_gravity_rotation_correction = s% use_gravity_rotation_correction
 eps_grav_factor = s% eps_grav_factor
 eps_mdot_factor = s% eps_mdot_factor
 include_composition_in_eps_grav = s% include_composition_in_eps_grav
 max_abs_rel_change_surf_lnS = s% max_abs_rel_change_surf_lnS
 max_num_surf_revisions = s% max_num_surf_revisions
 Gamma_lnS_eps_grav_full_off = s% Gamma_lnS_eps_grav_full_off
 Gamma_lnS_eps_grav_full_on = s% Gamma_lnS_eps_grav_full_on

 use_dPrad_dm_form_of_T_gradient_eqn = s% use_dPrad_dm_form_of_T_gradient_eqn
 use_gradT_actual_vs_gradT_MLT_for_T_gradient_eqn = s% use_gradT_actual_vs_gradT_MLT_for_T_gradient_eqn
 steps_before_use_velocity_time_centering = s% steps_before_use_velocity_time_centering
 include_P_in_velocity_time_centering = s% include_P_in_velocity_time_centering
 include_L_in_velocity_time_centering = s% include_L_in_velocity_time_centering
 P_theta_for_velocity_time_centering = s% P_theta_for_velocity_time_centering
 L_theta_for_velocity_time_centering = s% L_theta_for_velocity_time_centering
 use_P_d_1_div_rho_form_of_work_when_time_centering_velocity = s% use_P_d_1_div_rho_form_of_work_when_time_centering_velocity

 RTI_A = s% RTI_A
 RTI_B = s% RTI_B
 RTI_C = s% RTI_C
 RTI_D = s% RTI_D
 RTI_max_alpha = s% RTI_max_alpha
 RTI_C_X_factor = s% RTI_C_X_factor
 RTI_C_X0_frac = s% RTI_C_X0_frac
 RTI_dm_for_center_eta_nondecreasing = s% RTI_dm_for_center_eta_nondecreasing
 RTI_min_dm_behind_shock_for_full_on = s% RTI_min_dm_behind_shock_for_full_on
 RTI_energy_floor = s% RTI_energy_floor
 RTI_D_mix_floor = s% RTI_D_mix_floor
 RTI_min_m_for_D_mix_floor = s% RTI_min_m_for_D_mix_floor
 RTI_log_max_boost = s% RTI_log_max_boost 
 RTI_m_full_boost = s% RTI_m_full_boost
 RTI_m_no_boost = s% RTI_m_no_boost

 conv_vel_D = s% conv_vel_D
 conv_vel_siglimit = s% conv_vel_siglimit
 conv_vel_v0 = s% conv_vel_v0
 min_q_for_normal_mlt_gradT_full_off = s% min_q_for_normal_mlt_gradT_full_off
 max_q_for_normal_mlt_gradT_full_on = s% max_q_for_normal_mlt_gradT_full_on
 conv_vel_ignore_thermohaline = s% conv_vel_ignore_thermohaline
 conv_vel_ignore_semiconvection = s% conv_vel_ignore_semiconvection
 conv_vel_fully_lagrangian = s% conv_vel_fully_lagrangian
 conv_vel_include_homologous_term = s% conv_vel_include_homologous_term
 conv_vel_use_mlt_vc_start = s% conv_vel_use_mlt_vc_start

 velocity_logT_lower_bound = s% velocity_logT_lower_bound
 max_dt_yrs_for_velocity_logT_lower_bound = s% max_dt_yrs_for_velocity_logT_lower_bound
 velocity_q_upper_bound = s% velocity_q_upper_bound

 ! solvers

 tol_correction_norm = s% tol_correction_norm
 tol_max_correction = s% tol_max_correction
 correction_xa_limit = s% correction_xa_limit

 tol_correction_high_T_limit = s% tol_correction_high_T_limit
 tol_correction_norm_high_T = s% tol_correction_norm_high_T
 tol_max_correction_high_T = s% tol_max_correction_high_T

 tol_correction_extreme_T_limit = s% tol_correction_extreme_T_limit
 tol_correction_norm_extreme_T = s% tol_correction_norm_extreme_T
 tol_max_correction_extreme_T = s% tol_max_correction_extreme_T
 
 tol_bad_max_correction = s% tol_bad_max_correction
 bad_max_correction_series_limit = s% bad_max_correction_series_limit

 tol_residual_norm1 = s% tol_residual_norm1
 tol_max_residual1 = s% tol_max_residual1
 tol_residual_norm2 = s% tol_residual_norm2
 tol_max_residual2 = s% tol_max_residual2
 tol_residual_norm3 = s% tol_residual_norm3
 tol_max_residual3 = s% tol_max_residual3
 warning_limit_for_max_residual = s% warning_limit_for_max_residual
 trace_solver_damping = s% trace_solver_damping
 
 relax_use_gold_tolerances = s% relax_use_gold_tolerances
 relax_tol_correction_norm = s% relax_tol_correction_norm
 relax_tol_max_correction = s% relax_tol_max_correction
 relax_solver_iters_timestep_limit = s% relax_solver_iters_timestep_limit
 relax_iter_for_resid_tol2 = s% relax_iter_for_resid_tol2
 relax_tol_residual_norm1 = s% relax_tol_residual_norm1
 relax_tol_max_residual1 = s% relax_tol_max_residual1
 relax_iter_for_resid_tol3 = s% relax_iter_for_resid_tol3
 relax_tol_residual_norm2 = s% relax_tol_residual_norm2
 relax_tol_max_residual2 = s% relax_tol_max_residual2
 relax_tol_residual_norm3 = s% relax_tol_residual_norm3
 relax_tol_max_residual3 = s% relax_tol_max_residual3
 relax_maxT_for_gold_tolerances = s% relax_maxT_for_gold_tolerances
 
 use_gold_tolerances = s% use_gold_tolerances
 gold_solver_iters_timestep_limit = s% gold_solver_iters_timestep_limit 
 maxT_for_gold_tolerances = s% maxT_for_gold_tolerances
 gold_tol_residual_norm1 = s% gold_tol_residual_norm1
 gold_tol_max_residual1 = s% gold_tol_max_residual1
 gold_iter_for_resid_tol2 = s% gold_iter_for_resid_tol2
 gold_tol_residual_norm2 = s% gold_tol_residual_norm2
 gold_tol_max_residual2 = s% gold_tol_max_residual2
 gold_iter_for_resid_tol3 = s% gold_iter_for_resid_tol3
 gold_tol_residual_norm3 = s% gold_tol_residual_norm3
 gold_tol_max_residual3 = s% gold_tol_max_residual3
 steps_before_use_gold_tolerances = s% steps_before_use_gold_tolerances
 
 use_gold2_tolerances = s% use_gold2_tolerances
 gold2_solver_iters_timestep_limit = s% gold2_solver_iters_timestep_limit 
 gold2_tol_residual_norm1 = s% gold2_tol_residual_norm1
 gold2_tol_max_residual1 = s% gold2_tol_max_residual1
 gold2_iter_for_resid_tol2 = s% gold2_iter_for_resid_tol2
 gold2_tol_residual_norm2 = s% gold2_tol_residual_norm2
 gold2_tol_max_residual2 = s% gold2_tol_max_residual2
 gold2_iter_for_resid_tol3 = s% gold2_iter_for_resid_tol3
 gold2_tol_residual_norm3 = s% gold2_tol_residual_norm3
 gold2_tol_max_residual3 = s% gold2_tol_max_residual3
 steps_before_use_gold2_tolerances = s% steps_before_use_gold2_tolerances
 
 include_rotation_in_total_energy = s% include_rotation_in_total_energy

 convergence_ignore_equL_residuals = s% convergence_ignore_equL_residuals
 convergence_ignore_alpha_RTI_residuals = s% convergence_ignore_alpha_RTI_residuals

 iter_for_resid_tol2 = s% iter_for_resid_tol2
 iter_for_resid_tol3 = s% iter_for_resid_tol3

 solver_itermin = s% solver_itermin
 solver_itermin_until_reduce_min_corr_coeff = s% solver_itermin_until_reduce_min_corr_coeff
 solver_reduced_min_corr_coeff = s% solver_reduced_min_corr_coeff
 do_solver_damping_for_neg_xa = s% do_solver_damping_for_neg_xa
 hydro_mtx_max_allowed_abs_dlogT = s% hydro_mtx_max_allowed_abs_dlogT
 hydro_mtx_max_allowed_abs_dlogRho = s% hydro_mtx_max_allowed_abs_dlogRho
 min_logT_for_hydro_mtx_max_allowed = s% min_logT_for_hydro_mtx_max_allowed
 hydro_mtx_max_allowed_logT = s% hydro_mtx_max_allowed_logT
 hydro_mtx_max_allowed_logRho = s% hydro_mtx_max_allowed_logRho
 hydro_mtx_min_allowed_logT = s% hydro_mtx_min_allowed_logT
 hydro_mtx_min_allowed_logRho = s% hydro_mtx_min_allowed_logRho
 
 use_DGESVX_in_bcyclic = s% use_DGESVX_in_bcyclic
 use_equilibration_in_DGESVX = s% use_equilibration_in_DGESVX
 report_min_rcond_from_DGESXV = s% report_min_rcond_from_DGESXV
 
 op_split_burn = s% op_split_burn
 op_split_burn_min_T = s% op_split_burn_min_T
 op_split_burn_eps = s% op_split_burn_eps
 op_split_burn_odescal = s% op_split_burn_odescal
 op_split_burn_min_T_for_variable_T_solver = s% op_split_burn_min_T_for_variable_T_solver

 tiny_corr_coeff_limit = s% tiny_corr_coeff_limit
 scale_correction_norm = s% scale_correction_norm
 num_times_solver_reuse_mtx = s% num_times_solver_reuse_mtx
 corr_param_factor = s% corr_param_factor
 scale_max_correction = s% scale_max_correction
 ignore_min_corr_coeff_for_scale_max_correction = s% ignore_min_corr_coeff_for_scale_max_correction
 ignore_too_large_correction = s% ignore_too_large_correction
 ignore_species_in_max_correction = s% ignore_species_in_max_correction

 corr_norm_jump_limit = s% corr_norm_jump_limit
 max_corr_jump_limit = s% max_corr_jump_limit
 resid_norm_jump_limit = s% resid_norm_jump_limit
 max_resid_jump_limit = s% max_resid_jump_limit

 corr_coeff_limit = s% corr_coeff_limit
 tiny_corr_factor = s% tiny_corr_factor

 solver_max_tries_before_reject = s% solver_max_tries_before_reject
 max_tries1 = s% max_tries1
 max_tries_for_retry = s% max_tries_for_retry
 max_tries_after_5_retries = s% max_tries_after_5_retries
 max_tries_after_10_retries = s% max_tries_after_10_retries
 max_tries_after_20_retries = s% max_tries_after_20_retries
 retry_limit = s% retry_limit
 redo_limit = s% redo_limit

 use_Pvsc_art_visc = s% use_Pvsc_art_visc
 Pvsc_cq = s% Pvsc_cq
 Pvsc_zsh = s% Pvsc_zsh

 min_xa_hard_limit = s% min_xa_hard_limit
 min_xa_hard_limit_for_highT = s% min_xa_hard_limit_for_highT
 logT_max_for_min_xa_hard_limit = s% logT_max_for_min_xa_hard_limit
 logT_min_for_min_xa_hard_limit_for_highT = s% logT_min_for_min_xa_hard_limit_for_highT

 sum_xa_hard_limit = s% sum_xa_hard_limit
 sum_xa_hard_limit_for_highT = s% sum_xa_hard_limit_for_highT
 logT_max_for_sum_xa_hard_limit = s% logT_max_for_sum_xa_hard_limit
 logT_min_for_sum_xa_hard_limit_for_highT = s% logT_min_for_sum_xa_hard_limit_for_highT

 xa_clip_limit = s% xa_clip_limit
 report_solver_progress = s% report_solver_progress
 solver_test_partials_call_number = s% solver_test_partials_call_number
 solver_test_partials_iter_number = s% solver_test_partials_iter_number
 solver_epsder_chem = s% solver_epsder_chem
 solver_epsder_struct = s% solver_epsder_struct
 solver_numerical_jacobian = s% solver_numerical_jacobian
 solver_jacobian_nzlo = s% solver_jacobian_nzlo
 solver_jacobian_nzhi = s% solver_jacobian_nzhi
 solver_check_everything = s% solver_check_everything
 energy_conservation_dump_model_number = s% energy_conservation_dump_model_number
 solver_inspect_soln_flag = s% solver_inspect_soln_flag
 solver_test_partials_dx_0 = s% solver_test_partials_dx_0
 solver_test_partials_k = s% solver_test_partials_k
 solver_test_partials_k_low = s% solver_test_partials_k_low
 solver_test_partials_k_high = s% solver_test_partials_k_high
 solver_show_correction_info = s% solver_show_correction_info
 solver_test_partials_write_eos_call_info = s% solver_test_partials_write_eos_call_info
 solver_test_eos_partials = s% solver_test_eos_partials
 solver_test_kap_partials = s% solver_test_kap_partials
 solver_test_net_partials = s% solver_test_net_partials
 solver_test_atm_partials = s% solver_test_atm_partials
solver_test_partials_var_name = s% solver_test_partials_var_name
solver_test_partials_sink_name = s% solver_test_partials_sink_name
 solver_test_partials_equ_name = s% solver_test_partials_equ_name
 solver_test_partials_show_dx_var_name = s% solver_test_partials_show_dx_var_name
 solver_save_photo_call_number = s% solver_save_photo_call_number
 fill_arrays_with_NaNs = s% fill_arrays_with_NaNs
 zero_when_allocate = s% zero_when_allocate
 warn_when_large_rel_run_E_err = s% warn_when_large_rel_run_E_err
 warn_when_large_virial_thm_rel_err = s% warn_when_large_virial_thm_rel_err
 warn_when_get_a_bad_eos_result = s% warn_when_get_a_bad_eos_result
 warn_rates_for_high_temp = s% warn_rates_for_high_temp
 max_safe_logT_for_rates = s% max_safe_logT_for_rates
 eps_mdot_leak_frac_factor = s% eps_mdot_leak_frac_factor

 max_dt_div_tau_conv_for_TDC = s% max_dt_div_tau_conv_for_TDC
 max_dt_years_for_TDC = s% max_dt_years_for_TDC
 alpha_TDC_DAMP = s% alpha_TDC_DAMP
 alpha_TDC_DAMPR = s% alpha_TDC_DAMPR
 alpha_TDC_PtdVdt = s% alpha_TDC_PtdVdt
 max_X_for_gradT_eqn = s% max_X_for_gradT_eqn
 compare_TDC_to_MLT = s% compare_TDC_to_MLT

 RSP2_alfap= s% RSP2_alfap
 RSP2_alfad = s% RSP2_alfad
 RSP2_alfat= s% RSP2_alfat 
 RSP2_alfam= s% RSP2_alfam
 RSP2_alfar= s% RSP2_alfar
 RSP2_min_Lt_div_L_for_overshooting_mixing_type = s% RSP2_min_Lt_div_L_for_overshooting_mixing_type
 RSP2_min_Lc_div_L_for_convective_mixing_type = s% RSP2_min_Lc_div_L_for_convective_mixing_type
 RSP2_Lsurf_factor= s% RSP2_Lsurf_factor
 RSP2_use_Stellingwerf_Lr = s% RSP2_use_Stellingwerf_Lr
 RSP2_use_L_eqn_at_surface = s% RSP2_use_L_eqn_at_surface
 RSP2_assume_HSE = s% RSP2_assume_HSE
 RSP2_use_RSP_eqn_for_Y_face = s% RSP2_use_RSP_eqn_for_Y_face
 RSP2_use_mass_interp_face_values = s% RSP2_use_mass_interp_face_values
 RSP2_num_outermost_cells_forced_nonturbulent = s% RSP2_num_outermost_cells_forced_nonturbulent
 RSP2_num_innermost_cells_forced_nonturbulent = s% RSP2_num_innermost_cells_forced_nonturbulent
 RSP2_min_dt_div_tau_conv_switch_to_MLT = s% RSP2_min_dt_div_tau_conv_switch_to_MLT
 RSP2_min_dt_years_switch_to_MLT = s% RSP2_min_dt_years_switch_to_MLT
 RSP2_target_steps_per_cycle = s% RSP2_target_steps_per_cycle
 RSP2_max_num_periods = s% RSP2_max_num_periods
 RSP2_work_period = s% RSP2_work_period
 RSP2_map_first_period = s% RSP2_map_first_period
 RSP2_map_last_period = s% RSP2_map_last_period
 RSP2_min_max_R_for_periods = s% RSP2_min_max_R_for_periods
 RSP2_GREKM_avg_abs_frac_new = s% RSP2_GREKM_avg_abs_frac_new
 RSP2_GREKM_avg_abs_limit = s% RSP2_GREKM_avg_abs_limit
 RSP2_map_zone_interval = s% RSP2_map_zone_interval
 RSP2_work_filename = s% RSP2_work_filename
 RSP2_map_columns_filename = s% RSP2_map_columns_filename
 RSP2_map_filename = s% RSP2_map_filename
 RSP2_map_history_filename = s% RSP2_map_history_filename
 RSP2_write_map = s% RSP2_write_map
 RSP2_w_min_for_damping = s% RSP2_w_min_for_damping
 RSP2_source_seed = s% RSP2_source_seed
 RSP2_w_fix_if_neg = s% RSP2_w_fix_if_neg

 max_X_for_conv_timescale = s% max_X_for_conv_timescale
 min_X_for_conv_timescale = s% min_X_for_conv_timescale
 max_q_for_conv_timescale = s% max_q_for_conv_timescale
 min_q_for_conv_timescale = s% min_q_for_conv_timescale
 max_q_for_QHSE_timescale = s% max_q_for_QHSE_timescale
 min_q_for_QHSE_timescale = s% min_q_for_QHSE_timescale

 ! timestep
 max_timestep = s% max_timestep
 max_years_for_timestep = s% max_years_for_timestep

 hi_T_max_years_for_timestep = s% hi_T_max_years_for_timestep
 max_timestep_hi_T_limit = s% max_timestep_hi_T_limit

 min_timestep_factor = s% min_timestep_factor
 max_timestep_factor = s% max_timestep_factor
 max_timestep_factor_at_high_T = s% max_timestep_factor_at_high_T
 min_logT_for_max_timestep_factor_at_high_T = s% min_logT_for_max_timestep_factor_at_high_T
 time_delta_coeff = s% time_delta_coeff
 timestep_factor_for_retries = s% timestep_factor_for_retries
 retry_hold = s% retry_hold
 neg_mass_fraction_hold = s% neg_mass_fraction_hold
 timestep_dt_factor = s% timestep_dt_factor
 use_dt_low_pass_controller = s% use_dt_low_pass_controller
 
 force_timestep_min = s% force_timestep_min
 force_timestep_min_years = s% force_timestep_min_years
 force_timestep_min_factor = s% force_timestep_min_factor
 force_timestep = s% force_timestep
 force_timestep_years = s% force_timestep_years

 varcontrol_target = s% varcontrol_target
 min_allowed_varcontrol_target = s% min_allowed_varcontrol_target
 varcontrol_dt_limit_ratio_hard_max = s% varcontrol_dt_limit_ratio_hard_max
 xa_scale = s% xa_scale

 solver_iters_timestep_limit = s% solver_iters_timestep_limit

 burn_steps_limit = s% burn_steps_limit
 burn_steps_hard_limit = s% burn_steps_hard_limit

 diffusion_steps_limit = s% diffusion_steps_limit
 diffusion_steps_hard_limit = s% diffusion_steps_hard_limit
 diffusion_iters_limit = s% diffusion_iters_limit
 diffusion_iters_hard_limit = s% diffusion_iters_hard_limit

 dt_div_dt_cell_collapse_limit = s% dt_div_dt_cell_collapse_limit
 dt_div_dt_cell_collapse_hard_limit = s% dt_div_dt_cell_collapse_hard_limit
 dt_div_min_dr_div_cs_limit = s% dt_div_min_dr_div_cs_limit
 dt_div_min_dr_div_cs_hard_limit = s% dt_div_min_dr_div_cs_hard_limit
 
 min_abs_du_div_cs_for_dt_div_min_dr_div_cs_limit = s% min_abs_du_div_cs_for_dt_div_min_dr_div_cs_limit
 min_abs_u_div_cs_for_dt_div_min_dr_div_cs_limit = s% min_abs_u_div_cs_for_dt_div_min_dr_div_cs_limit
 min_k_for_dt_div_min_dr_div_cs_limit = s% min_k_for_dt_div_min_dr_div_cs_limit
 min_q_for_dt_div_min_dr_div_cs_limit = s% min_q_for_dt_div_min_dr_div_cs_limit
 max_q_for_dt_div_min_dr_div_cs_limit = s% max_q_for_dt_div_min_dr_div_cs_limit
 check_remnant_only_for_dt_div_min_dr_div_cs_limit = s% check_remnant_only_for_dt_div_min_dr_div_cs_limit

 dX_mix_dist_limit = s% dX_mix_dist_limit

 dH_limit_min_H = s% dH_limit_min_H
 dH_limit = s% dH_limit
 dH_hard_limit = s% dH_hard_limit
 dH_div_H_limit_min_H = s% dH_div_H_limit_min_H
 dH_div_H_limit = s% dH_div_H_limit
 dH_div_H_hard_limit = s% dH_div_H_hard_limit
 dH_decreases_only = s% dH_decreases_only

 dHe_limit_min_He = s% dHe_limit_min_He
 dHe_limit = s% dHe_limit
 dHe_hard_limit = s% dHe_hard_limit
 dHe_div_He_limit_min_He = s% dHe_div_He_limit_min_He
 dHe_div_He_limit = s% dHe_div_He_limit
 dHe_div_He_hard_limit = s% dHe_div_He_hard_limit
 dHe_decreases_only = s% dHe_decreases_only

 dHe3_limit_min_He3 = s% dHe3_limit_min_He3
 dHe3_limit = s% dHe3_limit
 dHe3_hard_limit = s% dHe3_hard_limit
 dHe3_div_He3_limit_min_He3 = s% dHe3_div_He3_limit_min_He3
 dHe3_div_He3_limit = s% dHe3_div_He3_limit
 dHe3_div_He3_hard_limit = s% dHe3_div_He3_hard_limit
 dHe3_decreases_only = s% dHe3_decreases_only

 dX_limit_min_X = s% dX_limit_min_X
 dX_limit = s% dX_limit
 dX_hard_limit = s% dX_hard_limit
 dX_div_X_limit_min_X = s% dX_div_X_limit_min_X
 dX_div_X_limit = s% dX_div_X_limit
 dX_div_X_hard_limit = s% dX_div_X_hard_limit
 dX_div_X_at_high_T_limit = s% dX_div_X_at_high_T_limit
 dX_div_X_at_high_T_hard_limit = s% dX_div_X_at_high_T_hard_limit
 dX_div_X_at_high_T_limit_lgT_min = s% dX_div_X_at_high_T_limit_lgT_min
 dX_decreases_only = s% dX_decreases_only

 dX_nuc_drop_min_X_limit = s% dX_nuc_drop_min_X_limit
 dX_nuc_drop_max_A_limit = s% dX_nuc_drop_max_A_limit
 dX_nuc_drop_limit = s% dX_nuc_drop_limit
 dX_nuc_drop_limit_at_high_T = s% dX_nuc_drop_limit_at_high_T
 dX_nuc_drop_hard_limit = s% dX_nuc_drop_hard_limit
 dX_nuc_drop_min_yrs_for_dt = s% dX_nuc_drop_min_yrs_for_dt

 dL_div_L_limit_min_L = s% dL_div_L_limit_min_L
 dL_div_L_limit = s% dL_div_L_limit
 dL_div_L_hard_limit = s% dL_div_L_hard_limit

 delta_lgP_limit = s% delta_lgP_limit
 delta_lgP_hard_limit = s% delta_lgP_hard_limit
 delta_lgP_limit_min_lgP = s% delta_lgP_limit_min_lgP

 delta_lgRho_limit = s% delta_lgRho_limit
 delta_lgRho_hard_limit = s% delta_lgRho_hard_limit
 delta_lgRho_limit_min_lgRho = s% delta_lgRho_limit_min_lgRho

 delta_lgT_limit = s% delta_lgT_limit
 delta_lgT_hard_limit = s% delta_lgT_hard_limit
 delta_lgT_limit_min_lgT = s% delta_lgT_limit_min_lgT

 delta_lgE_limit = s% delta_lgE_limit
 delta_lgE_hard_limit = s% delta_lgE_hard_limit
 delta_lgE_limit_min_lgE = s% delta_lgE_limit_min_lgE

 delta_lgR_limit = s% delta_lgR_limit
 delta_lgR_hard_limit = s% delta_lgR_hard_limit
 delta_lgR_limit_min_lgR = s% delta_lgR_limit_min_lgR

 delta_Ye_highT_limit = s% delta_Ye_highT_limit
 delta_Ye_highT_hard_limit = s% delta_Ye_highT_hard_limit
 minT_for_highT_Ye_limit = s% minT_for_highT_Ye_limit

 delta_lgL_nuc_cat_limit = s% delta_lgL_nuc_cat_limit
 delta_lgL_nuc_cat_hard_limit = s% delta_lgL_nuc_cat_hard_limit
 lgL_nuc_cat_burn_min = s% lgL_nuc_cat_burn_min
 lgL_nuc_mix_dist_limit = s% lgL_nuc_mix_dist_limit

 delta_lgL_H_limit = s% delta_lgL_H_limit
 delta_lgL_H_hard_limit = s% delta_lgL_H_hard_limit
 lgL_H_burn_min = s% lgL_H_burn_min
 lgL_H_drop_factor = s% lgL_H_drop_factor
 lgL_H_burn_relative_limit = s% lgL_H_burn_relative_limit

 delta_lgL_He_limit = s% delta_lgL_He_limit
 delta_lgL_He_hard_limit = s% delta_lgL_He_hard_limit
 lgL_He_burn_min = s% lgL_He_burn_min
 lgL_He_drop_factor = s% lgL_He_drop_factor
 lgL_He_burn_relative_limit = s% lgL_He_burn_relative_limit

 delta_lgL_z_limit = s% delta_lgL_z_limit
 delta_lgL_z_hard_limit = s% delta_lgL_z_hard_limit
 lgL_z_burn_min = s% lgL_z_burn_min
 lgL_z_drop_factor = s% lgL_z_drop_factor
 lgL_z_burn_relative_limit = s% lgL_z_burn_relative_limit

 delta_lgL_power_photo_limit = s% delta_lgL_power_photo_limit
 delta_lgL_power_photo_hard_limit = s% delta_lgL_power_photo_hard_limit
 lgL_power_photo_burn_min = s% lgL_power_photo_burn_min
 lgL_power_photo_drop_factor = s% lgL_power_photo_drop_factor
 min_lgT_for_lgL_power_photo_limit = s% min_lgT_for_lgL_power_photo_limit

 delta_lgL_nuc_limit = s% delta_lgL_nuc_limit
 delta_lgL_nuc_hard_limit = s% delta_lgL_nuc_hard_limit
 delta_lgL_nuc_at_high_T_limit = s% delta_lgL_nuc_at_high_T_limit
 delta_lgL_nuc_at_high_T_hard_limit = s% delta_lgL_nuc_at_high_T_hard_limit
 delta_lgL_nuc_at_high_T_limit_lgT_min = s% delta_lgL_nuc_at_high_T_limit_lgT_min
 
 max_lgT_for_lgL_nuc_limit = s% max_lgT_for_lgL_nuc_limit
 lgL_nuc_burn_min = s% lgL_nuc_burn_min
 lgL_nuc_drop_factor = s% lgL_nuc_drop_factor

 delta_lgRho_cntr_limit = s% delta_lgRho_cntr_limit
 delta_lgRho_cntr_hard_limit = s% delta_lgRho_cntr_hard_limit

 delta_lgT_cntr_limit = s% delta_lgT_cntr_limit
 delta_lgT_cntr_hard_limit = s% delta_lgT_cntr_hard_limit
 delta_lgT_cntr_limit_only_after_near_zams = s% delta_lgT_cntr_limit_only_after_near_zams

 delta_lgP_cntr_limit = s% delta_lgP_cntr_limit
 delta_lgP_cntr_hard_limit = s% delta_lgP_cntr_hard_limit

 delta_lgT_max_limit = s% delta_lgT_max_limit
 delta_lgT_max_hard_limit = s% delta_lgT_max_hard_limit
 delta_lgT_max_limit_lgT_min = s% delta_lgT_max_limit_lgT_min
 delta_lgT_max_limit_only_after_near_zams = s% delta_lgT_max_limit_only_after_near_zams

 delta_lgT_max_at_high_T_limit = s% delta_lgT_max_at_high_T_limit
 delta_lgT_max_at_high_T_hard_limit = s% delta_lgT_max_at_high_T_hard_limit
 delta_lgT_max_at_high_T_limit_lgT_min = s% delta_lgT_max_at_high_T_limit_lgT_min

 delta_log_eps_nuc_limit = s% delta_log_eps_nuc_limit
 delta_log_eps_nuc_hard_limit = s% delta_log_eps_nuc_hard_limit

 delta_dX_div_X_cntr_min = s% delta_dX_div_X_cntr_min
 delta_dX_div_X_cntr_max = s% delta_dX_div_X_cntr_max
 delta_dX_div_X_cntr_limit = s% delta_dX_div_X_cntr_limit
 delta_dX_div_X_cntr_hard_limit = s% delta_dX_div_X_cntr_hard_limit

 delta_dX_div_X_drop_only = s% delta_dX_div_X_drop_only
 delta_lg_XH_drop_only = s% delta_lg_XH_drop_only
 delta_lg_XHe_drop_only = s% delta_lg_XHe_drop_only
 delta_lg_XC_drop_only = s% delta_lg_XC_drop_only
 delta_lg_XNe_drop_only = s% delta_lg_XNe_drop_only
 delta_lg_XO_drop_only = s% delta_lg_XO_drop_only
 delta_lg_XSi_drop_only = s% delta_lg_XSi_drop_only
 delta_XH_drop_only = s% delta_XH_drop_only
 delta_XHe_drop_only = s% delta_XHe_drop_only
 delta_XC_drop_only = s% delta_XC_drop_only
 delta_XNe_drop_only = s% delta_XNe_drop_only
 delta_XO_drop_only = s% delta_XO_drop_only
 delta_XSi_drop_only = s% delta_XSi_drop_only

 delta_lg_XH_cntr_min = s% delta_lg_XH_cntr_min
 delta_lg_XH_cntr_max = s% delta_lg_XH_cntr_max
 delta_lg_XH_cntr_limit = s% delta_lg_XH_cntr_limit
 delta_lg_XH_cntr_hard_limit = s% delta_lg_XH_cntr_hard_limit

 delta_lg_XHe_cntr_min = s% delta_lg_XHe_cntr_min
 delta_lg_XHe_cntr_max = s% delta_lg_XHe_cntr_max
 delta_lg_XHe_cntr_limit = s% delta_lg_XHe_cntr_limit
 delta_lg_XHe_cntr_hard_limit = s% delta_lg_XHe_cntr_hard_limit

 delta_lg_XC_cntr_min = s% delta_lg_XC_cntr_min
 delta_lg_XC_cntr_max = s% delta_lg_XC_cntr_max
 delta_lg_XC_cntr_limit = s% delta_lg_XC_cntr_limit
 delta_lg_XC_cntr_hard_limit = s% delta_lg_XC_cntr_hard_limit

 delta_lg_XNe_cntr_limit = s% delta_lg_XNe_cntr_limit
 delta_lg_XNe_cntr_hard_limit = s% delta_lg_XNe_cntr_hard_limit
 delta_lg_XNe_cntr_min = s% delta_lg_XNe_cntr_min
 delta_lg_XNe_cntr_max = s% delta_lg_XNe_cntr_max

 delta_lg_XO_cntr_limit = s% delta_lg_XO_cntr_limit
 delta_lg_XO_cntr_hard_limit = s% delta_lg_XO_cntr_hard_limit
 delta_lg_XO_cntr_min = s% delta_lg_XO_cntr_min
 delta_lg_XO_cntr_max = s% delta_lg_XO_cntr_max

 delta_lg_XSi_cntr_limit = s% delta_lg_XSi_cntr_limit
 delta_lg_XSi_cntr_hard_limit = s% delta_lg_XSi_cntr_hard_limit
 delta_lg_XSi_cntr_min = s% delta_lg_XSi_cntr_min
 delta_lg_XSi_cntr_max = s% delta_lg_XSi_cntr_max

 delta_XH_cntr_limit = s% delta_XH_cntr_limit
 delta_XH_cntr_hard_limit = s% delta_XH_cntr_hard_limit
 delta_XHe_cntr_limit = s% delta_XHe_cntr_limit
 delta_XHe_cntr_hard_limit = s% delta_XHe_cntr_hard_limit
 delta_XC_cntr_limit = s% delta_XC_cntr_limit
 delta_XC_cntr_hard_limit = s% delta_XC_cntr_hard_limit
 delta_XNe_cntr_limit = s% delta_XNe_cntr_limit
 delta_XNe_cntr_hard_limit = s% delta_XNe_cntr_hard_limit
 delta_XO_cntr_limit = s% delta_XO_cntr_limit
 delta_XO_cntr_hard_limit = s% delta_XO_cntr_hard_limit
 delta_XSi_cntr_limit = s% delta_XSi_cntr_limit
 delta_XSi_cntr_hard_limit = s% delta_XSi_cntr_hard_limit

 delta_lgTeff_limit = s% delta_lgTeff_limit
 delta_lgTeff_hard_limit = s% delta_lgTeff_hard_limit

 delta_lgL_limit = s% delta_lgL_limit
 delta_lgL_limit_L_min = s% delta_lgL_limit_L_min
 delta_lgL_hard_limit = s% delta_lgL_hard_limit

 delta_HR_ds_L = s% delta_HR_ds_L
 delta_HR_ds_Teff = s% delta_HR_ds_Teff
 delta_HR_limit = s% delta_HR_limit
 delta_HR_hard_limit = s% delta_HR_hard_limit

 delta_lg_star_mass_limit = s% delta_lg_star_mass_limit
 delta_lg_star_mass_hard_limit = s% delta_lg_star_mass_hard_limit

 delta_mdot_atol = s% delta_mdot_atol
 delta_mdot_rtol = s% delta_mdot_rtol
 delta_mdot_limit = s% delta_mdot_limit
 delta_mdot_hard_limit = s% delta_mdot_hard_limit

 adjust_J_q_limit = s% adjust_J_q_limit
 adjust_J_q_hard_limit = s% adjust_J_q_hard_limit
 never_skip_hard_limits = s% never_skip_hard_limits
 relax_hard_limits_after_retry = s% relax_hard_limits_after_retry
 report_dt_hard_limit_retries = s% report_dt_hard_limit_retries
 report_min_dr_div_cs = s% report_min_dr_div_cs
 report_solver_dt_info = s% report_solver_dt_info

 limit_for_rel_error_in_energy_conservation = s% limit_for_rel_error_in_energy_conservation
 hard_limit_for_rel_error_in_energy_conservation = s% hard_limit_for_rel_error_in_energy_conservation

 min_chem_eqn_scale = s% min_chem_eqn_scale

 ! controls for the evolve routine
 trace_evolve = s% trace_evolve


 ! misc
 zams_filename = s% zams_filename
 set_rho_to_dm_div_dV = s% set_rho_to_dm_div_dV

 use_other_eos = s% use_other_eos
 use_other_surface_PT = s% use_other_surface_PT
 use_other_kap = s% use_other_kap
 use_other_diffusion = s% use_other_diffusion
 use_other_diffusion_factor = s% use_other_diffusion_factor
 use_other_adjust_mdot = s% use_other_adjust_mdot
 use_other_j_for_adjust_J_lost = s% use_other_j_for_adjust_J_lost
 use_other_alpha_mlt = s% use_other_alpha_mlt
 use_other_am_mixing = s% use_other_am_mixing
 use_other_brunt = s% use_other_brunt
 use_other_brunt_smoothing = s% use_other_brunt_smoothing
 use_other_solver_monitor = s% use_other_solver_monitor
 use_other_build_initial_model = s% use_other_build_initial_model
 use_other_cgrav = s% use_other_cgrav
 use_other_mesh_delta_coeff_factor = s% use_other_mesh_delta_coeff_factor
 use_other_energy_implicit = s% use_other_energy_implicit
 use_other_momentum_implicit = s% use_other_momentum_implicit
 use_other_momentum = s% use_other_momentum
 use_other_remove_surface = s% use_other_remove_surface
 use_other_energy = s% use_other_energy
 use_other_pressure = s% use_other_pressure
 use_other_mesh_functions = s% use_other_mesh_functions
 use_other_eps_grav = s% use_other_eps_grav
 use_other_gradr_factor = s% use_other_gradr_factor
 use_other_D_mix = s% use_other_D_mix
 use_other_neu = s% use_other_neu
 use_other_net_get = s% use_other_net_get
 use_other_opacity_factor = s% use_other_opacity_factor
 use_other_diffusion_coefficients = s% use_other_diffusion_coefficients
 use_other_pgstar_plots = s% use_other_pgstar_plots
 use_other_eval_fp_ft = s% use_other_eval_fp_ft
 use_other_eval_i_rot = s% use_other_eval_i_rot
 use_other_torque = s% use_other_torque
 use_other_torque_implicit = s% use_other_torque_implicit
 use_other_wind = s% use_other_wind
 use_other_accreting_state = s% use_other_accreting_state
 use_other_after_struct_burn_mix = s% use_other_after_struct_burn_mix
 use_other_before_struct_burn_mix = s% use_other_before_struct_burn_mix
 use_other_astero_freq_corr = s% use_other_astero_freq_corr
 use_other_timestep_limit = s% use_other_timestep_limit
 use_other_set_pgstar_controls = s% use_other_set_pgstar_controls
 use_other_screening = s% use_other_screening

 x_ctrl = s% x_ctrl
 x_integer_ctrl = s% x_integer_ctrl
 x_logical_ctrl = s% x_logical_ctrl
 x_character_ctrl = s% x_character_ctrl

 ! info for debugging
 stop_for_bad_nums = s% stop_for_bad_nums
 report_ierr = s% report_ierr
 report_bad_negative_xa = s% report_bad_negative_xa

 diffusion_dump_call_number = s% diffusion_dump_call_number


 end subroutine set_controls_for_writing

 end module ctrls_io

