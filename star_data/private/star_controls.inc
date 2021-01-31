


      ! starting specifications

         real(dp) :: initial_mass
         real(dp) :: initial_z
         real(dp) :: initial_y
         real(dp) :: initial_he3

      ! definition of core boundaries
         real(dp) :: &
            he_core_boundary_h1_fraction, &
            co_core_boundary_he4_fraction, &
            fe_core_boundary_si28_fraction, &
            neutron_rich_core_boundary_Ye_max, &
            min_boundary_fraction

      ! when to stop

         integer :: max_model_number
         real(dp) :: when_to_stop_atol
         real(dp) :: when_to_stop_rtol

         real(dp) :: max_age, max_age_in_seconds, &
            dt_years_for_steps_before_max_age, reduction_factor_for_max_timestep
         integer :: num_adjusted_dt_steps_before_max_age
         real(dp) :: gamma_center_limit
         real(dp) :: eta_center_limit
         real(dp) :: log_center_temp_limit
         real(dp) :: log_center_temp_lower_limit
         real(dp) :: log_center_density_limit, log_center_density_lower_limit
         real(dp) :: log_max_temp_upper_limit, log_max_temp_lower_limit
         real(dp) :: surface_accel_div_grav_limit
         real(dp) :: center_entropy_limit, center_entropy_lower_limit
         real(dp) :: max_entropy_limit, max_entropy_lower_limit

         character (len=iso_name_length) :: &
            xa_central_lower_limit_species(num_xa_central_limits)
         real(dp) :: xa_central_lower_limit(num_xa_central_limits)

         character (len=iso_name_length) :: &
            xa_central_upper_limit_species(num_xa_central_limits)
         real(dp) :: xa_central_upper_limit(num_xa_central_limits)

         character (len=iso_name_length) :: &
            xa_surface_lower_limit_species(num_xa_surface_limits)
         real(dp) :: xa_surface_lower_limit(num_xa_surface_limits)

         character (len=iso_name_length) :: &
            xa_surface_upper_limit_species(num_xa_surface_limits)
         real(dp) :: xa_surface_upper_limit(num_xa_surface_limits)

         character (len=iso_name_length) :: &
            xa_average_lower_limit_species(num_xa_average_limits)
         real(dp) :: xa_average_lower_limit(num_xa_average_limits)

         character (len=iso_name_length) :: &
            xa_average_upper_limit_species(num_xa_average_limits)
         real(dp) :: xa_average_upper_limit(num_xa_average_limits)

         real(dp) :: HB_limit

         real(dp) :: star_mass_min_limit, star_mass_max_limit
         real(dp) :: bound_mass_min_limit, bound_mass_max_limit
         
         real(dp) :: star_species_mass_min_limit, star_species_mass_max_limit
         character (len=iso_name_length) :: &
            star_species_mass_min_limit_iso, star_species_mass_max_limit_iso

         real(dp) :: xmstar_min_limit, xmstar_max_limit

         real(dp) :: envelope_mass_limit
         real(dp) :: envelope_fraction_left_limit

         real(dp) :: he_core_mass_limit
         real(dp) :: co_core_mass_limit
         real(dp) :: fe_core_mass_limit
         real(dp) :: neutron_rich_core_mass_limit

         real(dp) :: he_layer_mass_lower_limit

         real(dp) :: abs_diff_lg_LH_lg_Ls_limit

         real(dp) :: Teff_upper_limit
         real(dp) :: Teff_lower_limit

         real(dp) :: photosphere_r_upper_limit
         real(dp) :: photosphere_r_lower_limit

         real(dp) :: photosphere_m_upper_limit
         real(dp) :: photosphere_m_lower_limit
         real(dp) :: photosphere_m_sub_M_center_limit

         real(dp) :: log_Teff_upper_limit
         real(dp) :: log_Teff_lower_limit

         real(dp) :: log_Tsurf_upper_limit
         real(dp) :: log_Tsurf_lower_limit

         real(dp) :: log_Rsurf_upper_limit
         real(dp) :: log_Rsurf_lower_limit

         real(dp) :: log_L_upper_limit
         real(dp) :: log_L_lower_limit

         real(dp) :: log_g_upper_limit
         real(dp) :: log_g_lower_limit

         real(dp) :: log_Psurf_upper_limit
         real(dp) :: log_Psurf_lower_limit

         real(dp) :: log_Dsurf_upper_limit
         real(dp) :: log_Dsurf_lower_limit

         real(dp) :: power_nuc_burn_upper_limit
         real(dp) :: power_h_burn_upper_limit
         real(dp) :: power_he_burn_upper_limit
         real(dp) :: power_z_burn_upper_limit

         real(dp) :: power_nuc_burn_lower_limit
         real(dp) :: power_h_burn_lower_limit
         real(dp) :: power_he_burn_lower_limit
         real(dp) :: power_z_burn_lower_limit

         integer :: max_number_retries, relax_max_number_retries
         real(dp) :: max_abs_rel_run_E_err

         real(dp) :: min_timestep_limit
         real(dp) :: center_Ye_lower_limit, center_R_lower_limit, &
            fe_core_infall_limit, non_fe_core_infall_limit, non_fe_core_rebound_limit
         real(dp) :: v_div_csound_surf_limit, v_div_csound_max_limit
         real(dp) :: Lnuc_div_L_upper_limit, Lnuc_div_L_lower_limit
         real(dp) :: v_surf_div_v_kh_upper_limit, v_surf_div_v_kh_lower_limit
         real(dp) :: v_surf_div_v_esc_limit, Lnuc_div_L_zams_limit, v_surf_kms_limit
         logical :: stop_near_zams, stop_at_phase_PreMS, & 
            stop_at_phase_ZAMS, stop_at_phase_IAMS, stop_at_phase_TAMS, &
            stop_at_phase_He_Burn, stop_at_phase_ZACHeB, stop_at_phase_TACHeB, &
            stop_at_phase_TP_AGB, stop_at_phase_C_Burn, stop_at_phase_Ne_Burn, &
            stop_at_phase_O_Burn, stop_at_phase_Si_Burn, stop_at_phase_WDCS
         real(dp) :: Pgas_div_P_limit, Pgas_div_P_limit_max_q
         real(dp) :: gamma1_limit, gamma1_limit_max_q, gamma1_limit_max_v_div_vesc
         real(dp) :: peak_burn_vconv_div_cs_limit, omega_div_omega_crit_limit
         real(dp) :: delta_nu_lower_limit, delta_nu_upper_limit
         real(dp) :: delta_Pg_lower_limit, delta_Pg_upper_limit
         real(dp) :: shock_mass_upper_limit, mach1_mass_upper_limit, &
            stop_when_reach_this_cumulative_extra_heating

      ! output of "snapshots" for restarts

         integer :: photo_interval, photo_digits
         character (len=strlen) :: photo_directory


      ! output of logs and profiles

         logical :: do_history_file
         integer :: history_interval
         integer :: terminal_interval
         integer :: write_header_frequency
         character (len=strlen) :: extra_terminal_output_file, &
            terminal_show_age_units, terminal_show_timestep_units
         logical :: terminal_show_log_dt, terminal_show_log_age
         integer :: num_trace_history_values
         character (len=maxlen_history_column_name) :: trace_history_value_name(50)

         character (len=strlen) :: log_directory
         character (len=strlen) :: star_history_name
         character (len=strlen) :: star_history_header_name
         character (len=strlen) :: star_history_dbl_format, &
            star_history_int_format, star_history_txt_format

         logical :: write_profiles_flag, profile_header_include_sys_details
         integer :: profile_interval, priority_profile_interval
         character (len=strlen) :: profiles_index_name
         character (len=strlen) :: profile_data_prefix
         character (len=strlen) :: profile_data_suffix
         character (len=strlen) :: profile_data_header_suffix
         character (len=strlen) :: profile_dbl_format, profile_int_format, profile_txt_format

         integer :: profile_model
         integer :: max_num_profile_models, max_num_profile_zones

         logical :: write_controls_info_with_profile
         character (len=strlen) :: controls_data_prefix, controls_data_suffix

         logical :: write_model_with_profile
         character (len=strlen) :: model_data_prefix, model_data_suffix

         logical :: write_pulse_data_with_profile
         character (len=strlen) :: pulse_data_format, &
            format_for_FGONG_data, format_for_OSC_data
         logical :: &
            add_atmosphere_to_pulse_data, &
            add_center_point_to_pulse_data, &
            keep_surface_point_for_pulse_data, &
            add_double_points_to_pulse_data, &
            interpolate_rho_for_pulse_data
         real(dp) :: threshold_grad_mu_for_double_point
         integer :: max_number_of_double_points

         integer :: max_num_gyre_points

         real(dp) :: fgong_zero_A_inside_r
         real(dp) :: trace_mass_location, &
            min_tau_for_max_abs_v_location, min_q_for_inner_mach1_location, &
            max_q_for_outer_mach1_location, mass_depth_for_L_surf, &
            conv_core_gap_dq_limit
         real(dp) :: mixing_D_limit_for_log

      ! burn zone eps definitions for use in logs and profiles

         real(dp) :: burn_min1
         real(dp) :: burn_min2

      ! limit for max conv vel div csound

         real(dp) :: max_conv_vel_div_csound_maxq, max_q_for_limit_conv_vel, &
            max_mass_in_gm_for_limit_conv_vel, max_r_in_cm_for_limit_conv_vel
         integer :: width_for_limit_conv_vel

      ! for reported averages

         real(dp) :: surface_avg_abundance_dq
         real(dp) :: center_avg_value_dq

      ! mixing regions

         logical :: remove_mixing_glitches, okay_to_remove_mixing_singleton, &
            redo_conv_for_dr_lt_mixing_length
         real(dp) :: min_convective_gap
         real(dp) :: min_thermohaline_gap
         real(dp) :: min_semiconvection_gap
         real(dp) :: min_thermohaline_dropout
         real(dp) :: max_dropout_gradL_sub_grada
         logical :: remove_embedded_semiconvection
         logical  :: recalc_mix_info_after_evolve
         real(dp) :: prune_bad_cz_min_Hp_height
         real(dp) :: prune_bad_cz_min_log_eps_nuc


      ! mixing parameters

         real(dp) :: alpha_semiconvection, semiconvection_upper_limit_center_h1
         real(dp) :: thermohaline_coeff
         character (len=32) :: &
            thermohaline_option, MLT_option, semiconvection_option
         real(dp) :: mixing_length_alpha, remove_small_D_limit
         logical :: use_Ledoux_criterion
         integer :: num_cells_for_smooth_gradL_composition_term
         real(dp) :: threshold_for_smooth_gradL_composition_term
         real(dp) :: clip_D_limit
         logical :: alt_scale_height_flag, fix_eps_grav_transition_to_grid
         real(dp) :: Henyey_MLT_y_param, Henyey_MLT_nu_param, mlt_Pturb_factor, &
            min_logT_for_make_gradr_sticky_in_solver_iters
         logical :: mlt_use_rotation_correction, make_gradr_sticky_in_solver_iters, &
            no_MLT_below_shock, no_MLT_below_T_max

         real(dp) :: T_mix_limit

         logical :: okay_to_reduce_gradT_excess
         real(dp) :: &
            gradT_excess_f1, &
            gradT_excess_f2, &
            gradT_excess_max_center_h1, &
            gradT_excess_min_center_he4, &
            gradT_excess_max_logT, &
            gradT_excess_min_log_tau_full_on, &
            gradT_excess_max_log_tau_full_off

         real(dp) :: mlt_gradT_fraction, max_logT_for_mlt

      ! mixing parameters

         real(dp) :: min_overshoot_q, overshoot_alpha
         
         real(dp) :: burn_z_mix_region_logT, burn_he_mix_region_logT, burn_h_mix_region_logT
         real(dp) :: max_Y_for_burn_z_mix_region, max_X_for_burn_he_mix_region
         
         logical :: limit_overshoot_Hp_using_size_of_convection_zone

         logical        :: predictive_mix(NUM_PREDICTIVE_PARAM_SETS)
         real(dp)       :: predictive_superad_thresh(NUM_PREDICTIVE_PARAM_SETS)
         character(256) :: predictive_avoid_reversal(NUM_PREDICTIVE_PARAM_SETS)
         character(256) :: predictive_limit_ingestion(NUM_PREDICTIVE_PARAM_SETS)
         real(dp)       :: predictive_ingestion_factor(NUM_PREDICTIVE_PARAM_SETS)
         character(256) :: predictive_zone_type(NUM_PREDICTIVE_PARAM_SETS)
         character(256) :: predictive_zone_loc(NUM_PREDICTIVE_PARAM_SETS)
         character(256) :: predictive_bdy_loc(NUM_PREDICTIVE_PARAM_SETS)
         real(dp)       :: predictive_bdy_q_min(NUM_PREDICTIVE_PARAM_SETS)
         real(dp)       :: predictive_bdy_q_max(NUM_PREDICTIVE_PARAM_SETS)

         logical        :: do_conv_premix
         logical        :: conv_premix_avoid_increase
         real(dp)       :: conv_premix_time_factor
         logical        :: conv_premix_fix_pgas
         logical        :: conv_premix_dump_snapshots
         logical        :: do_premix_heating ! include eps_pre_mix in energy equation

         real(dp)       :: overshoot_f(NUM_OVERSHOOT_PARAM_SETS)
         real(dp)       :: overshoot_f0(NUM_OVERSHOOT_PARAM_SETS)
         real(dp)       :: overshoot_D0(NUM_OVERSHOOT_PARAM_SETS)
         real(dp)       :: overshoot_Delta0(NUM_OVERSHOOT_PARAM_SETS)
         real(dp)       :: overshoot_mass_full_on(NUM_OVERSHOOT_PARAM_SETS)
         real(dp)       :: overshoot_mass_full_off(NUM_OVERSHOOT_PARAM_SETS)
         character(256) :: overshoot_zone_type(NUM_OVERSHOOT_PARAM_SETS)
         character(256) :: overshoot_zone_loc(NUM_OVERSHOOT_PARAM_SETS)
         character(256) :: overshoot_bdy_loc(NUM_OVERSHOOT_PARAM_SETS)
         character(256) :: overshoot_scheme(NUM_OVERSHOOT_PARAM_SETS)
         real(dp)       :: overshoot_D_min
         real(dp)       :: overshoot_brunt_B_max   

         logical :: RSP_write_map, RSP_report_limit_dt, &
            use_other_RSP_linear_analysis, RSP_testing, RSP_use_Prad_for_Psurf, RSP_use_atm_grey_with_kap_for_Psurf, &
            RSP_fixed_Psurf, RSP_hydro_only, RSP_report_undercorrections, RSP_relax_initial_model, &
            use_other_RSP_build_model, RSP_relax_adjust_inner_mass_distribution, RSP_trace_RSP_build_model, &
            RSP_relax_alfap_before_alfat, use_RSP_new_start_scheme
         integer :: RSP_max_num_periods, RSP_max_iters_per_try, RSP_max_retries_per_step, &
            RSP_mode_for_setting_PERIODLIN, RSP_nmodes, &
            RSP_target_steps_per_cycle, RSP_nz_outer, RSP_nz, &
            RSP_map_first_period, RSP_map_last_period, RSP_map_zone_interval, RSP_relax_max_tries, &
            RSP_max_outer_dm_tries, RSP_work_period, RSP_max_inner_scale_tries
         character(256) :: RSP_map_filename, RSP_map_history_filename, &
            RSP_work_filename, RSP_map_columns_filename
         real(dp) :: &
            RSP_theta, RSP_thetat, RSP_thetaq, RSP_thetau, RSP_thetae, RSP_wtr, RSP_wtc, RSP_wtt, RSP_gam, &
            RSP_alfa, RSP_alfap, RSP_alfam, RSP_alfat, RSP_alfas, RSP_alfac, RSP_alfad, RSP_gammar, RSP_efl0, &
            RSP_cq, RSP_zsh, RSP_Qvisc_quadratic, RSP_Qvisc_linear, RSP_tol_max_corr, RSP_tol_max_resid, &
            RSP_nz_div_IBOTOM, RSP_Teff, RSP_mass, RSP_L, RSP_X, RSP_Z, RSP_GREKM_avg_abs_frac_new, &
            RSP_min_max_R_for_periods, RSP_GREKM_avg_abs_limit, RSP_Avel, RSP_Arnd, RSP_default_PERIODLIN, &
            RSP_kick_vsurf_km_per_sec, RSP_T_anchor, RSP_T_inner, RSP_Psurf, set_RSP_Psurf_to_multiple_of_initial_P1, &
            RSP_min_PERIOD_div_PERIODLIN, RSP_v_div_cs_threshold_for_dt_limit, RSP_initial_dt_factor, &
            RSP_fraction_1st_overtone, RSP_fraction_2nd_overtone, &
            RSP_max_dt_times_min_dr_div_cs, RSP_tau_surf_for_atm_grey_with_kap, &
            RSP_min_tau_for_turbulent_flux, RSP_max_dt, &
            RSP_kap_density_factor, RSP_T_anchor_tolerance, RSP_T_inner_tolerance, RSP_relax_dm_tolerance, &
            RSP_max_dt_times_min_rad_diff_time, RSP_Qvisc_linear_static, RSP_dq_1_factor, &
            RSP_surface_tau
         
         real(dp) :: &
            RSP_min_Y, &
            RSP_alpha, RSP_alpha_m, RSP_alpha_t, RSP_alpha_P_standard, &
            RSP_alpha_s_standard, RSP_alpha_c_standard, &
            RSP_alpha_d_standard, RSP_gamma_r_standard, &
            RSP_alpha_P_factor, RSP_alpha_s_factor, RSP_alpha_c_factor, &
            RSP_alpha_d_factor, RSP_gamma_r_factor, &
            RSP_min_tau, RSP_min_q, max_conv_vel_div_csound, &
            min_dt_for_increases_in_convection_velocity, &
            max_v_for_convection, max_q_for_convection_with_hydro_on, &
            max_v_div_cs_for_convection, &
            max_abs_du_div_cs_for_convection, &
            RTI_smooth_fraction, RTI_smooth_mass, &
            alpha_RTI_diffusion_factor, &
            dudt_RTI_diffusion_factor, &
            dedt_RTI_diffusion_factor, &
            dlnddt_RTI_diffusion_factor, &
            composition_RTI_diffusion_factor, &
            max_M_RTI_factors_full_on, min_M_RTI_factors_full_off, &
            alpha_RTI_src_min_v_div_cs, &
            alpha_RTI_src_max_q, alpha_RTI_src_min_q
         integer :: RTI_smooth_iterations

      ! brunt
         logical :: calculate_Brunt_B, calculate_Brunt_N2, use_brunt_gradmuX_form
         real(dp) :: min_magnitude_brunt_B, brunt_N2_coefficient
         integer :: num_cells_for_smooth_brunt_B
         real(dp) :: threshold_for_smooth_brunt_B

      ! atmosphere -- surface boundary conditions

         character (len=strlen) :: atm_option
         character (len=strlen) :: atm_off_table_option
         real(dp) :: Pextra_factor
         real(dp) :: atm_fixed_Teff
         real(dp) :: atm_fixed_Tsurf
         real(dp) :: atm_fixed_Psurf

         character (len=strlen) :: atm_T_tau_relation
         character (len=strlen) :: atm_T_tau_opacity
         real(dp) :: atm_T_tau_errtol
         integer :: atm_T_tau_max_iters
         integer :: atm_T_tau_max_steps

         character (len=strlen) :: atm_table

         character (len=strlen) :: atm_irradiated_opacity
         real(dp) :: atm_irradiated_errtol
         real(dp) :: atm_irradiated_T_eq
         real(dp) :: atm_irradiated_kap_v
         real(dp) :: atm_irradiated_kap_v_div_kap_th
         real(dp) :: atm_irradiated_P_surf
         integer :: atm_irradiated_max_iters

         real(dp) :: surface_extra_Pgas
         real(dp) :: Tsurf_factor
         real(dp) :: fixed_vsurf
         real(dp) :: fixed_L_for_BB_outer_BC
         real(dp) :: tau_for_L_BB
         logical :: use_atm_PT_at_center_of_surface_cell
         logical :: use_compression_outer_BC
         logical :: use_momentum_outer_BC
         logical :: use_zero_Pgas_outer_BC
         logical :: use_zero_dLdm_outer_BC
         logical :: use_T_black_body_outer_BC
         logical :: use_fixed_vsurf_outer_BC
         logical :: use_fixed_L_for_BB_outer_BC

      ! create atmosphere

         real(dp) :: atm_build_tau_outer
         real(dp) :: atm_build_dlogtau
         real(dp) :: atm_build_errtol

         logical :: use_T_tau_gradr_factor

      ! mass gain or loss

         real(dp) :: mass_change

         real(dp) :: mass_change_full_on_dt, mass_change_full_off_dt
         logical :: trace_dt_control_mass_change

         logical :: use_accreted_material_j
         real(dp) :: accreted_material_j

         real(dp) :: cool_wind_full_on_T
         real(dp) :: hot_wind_full_on_T

         real(dp) :: min_wind, max_wind

         real(dp) :: mdot_omega_power, max_rotational_mdot_boost, lim_trace_rotational_mdot_boost, &
            rotational_mdot_boost_fac, rotational_mdot_kh_fac, surf_avg_tau, surf_avg_tau_min, &
            max_mdot_jump_for_rotation
         
         integer :: max_tries_for_implicit_wind
         real(dp) :: iwind_tolerance, iwind_lambda

         real(dp) :: super_eddington_scaling_factor, super_eddington_wind_Ledd_factor, &
            wind_boost_full_off_L_div_Ledd, wind_boost_full_on_L_div_Ledd, &
            super_eddington_wind_max_boost
         logical :: trace_super_eddington_wind_boost

         character (len=strlen) :: hot_wind_scheme
         character (len=strlen) :: cool_wind_RGB_scheme
         character (len=strlen) :: cool_wind_AGB_scheme
         real(dp) :: RGB_to_AGB_wind_switch

         real(dp) :: Reimers_scaling_factor
         real(dp) :: Blocker_scaling_factor
         real(dp) :: de_Jager_scaling_factor
         real(dp) :: van_Loon_scaling_factor
         real(dp) :: Nieuwenhuijzen_scaling_factor
         real(dp) :: Vink_scaling_factor
         real(dp) :: Dutch_scaling_factor
         character (len=strlen) :: Dutch_wind_lowT_scheme

         real(dp) :: wind_H_envelope_limit, wind_H_He_envelope_limit, wind_He_layer_limit

         logical :: no_wind_if_no_rotation

         real(dp) :: rlo_scaling_factor, rlo_wind_min_L, rlo_wind_max_Teff
         real(dp) :: rlo_wind_roche_lobe_radius, rlo_wind_base_mdot
         real(dp) :: roche_lobe_xfer_full_on, roche_lobe_xfer_full_off
         real(dp) :: rlo_wind_scale_height

         real(dp) :: max_logT_for_k_below_const_q, max_logT_for_k_const_mass
         real(dp) :: min_q_for_k_below_const_q, max_q_for_k_below_const_q
         real(dp) :: min_q_for_k_const_mass, max_q_for_k_const_mass

         ! composition of accreted material
         logical :: accrete_same_as_surface

         logical :: accrete_given_mass_fractions
         integer :: num_accretion_species
         character (len=iso_name_length) :: accretion_species_id(max_num_accretion_species)
         real(dp) :: accretion_species_xa(max_num_accretion_species)

         real(dp) :: accretion_h1, accretion_h2, accretion_he3, accretion_he4
         integer :: accretion_zfracs
         logical :: accretion_dump_missing_metals_into_heaviest

         ! special list of z fractions
         real(dp) :: z_fraction_li
         real(dp) :: z_fraction_be
         real(dp) :: z_fraction_b
         real(dp) :: z_fraction_c
         real(dp) :: z_fraction_n
         real(dp) :: z_fraction_o
         real(dp) :: z_fraction_f
         real(dp) :: z_fraction_ne
         real(dp) :: z_fraction_na
         real(dp) :: z_fraction_mg
         real(dp) :: z_fraction_al
         real(dp) :: z_fraction_si
         real(dp) :: z_fraction_p
         real(dp) :: z_fraction_s
         real(dp) :: z_fraction_cl
         real(dp) :: z_fraction_ar
         real(dp) :: z_fraction_k
         real(dp) :: z_fraction_ca
         real(dp) :: z_fraction_sc
         real(dp) :: z_fraction_ti
         real(dp) :: z_fraction_v
         real(dp) :: z_fraction_cr
         real(dp) :: z_fraction_mn
         real(dp) :: z_fraction_fe
         real(dp) :: z_fraction_co
         real(dp) :: z_fraction_ni
         real(dp) :: z_fraction_cu
         real(dp) :: z_fraction_zn

         real(dp) :: lgT_hi_for_set_new_abundances, lgT_lo_for_set_new_abundances


      ! automatic stops for mass loss/gain

         real(dp) :: max_star_mass_for_gain, min_star_mass_for_loss
         real(dp) :: max_T_center_for_any_mass_loss, max_T_center_for_full_mass_loss


      ! extra power source

         real(dp) :: extra_power_source


      ! relaxation parameters

         real(dp) :: relax_dY
         real(dp) :: relax_dlnZ


      ! mesh adjustment controls

         real(dp) :: mesh_delta_coeff
         real(dp) :: mesh_delta_coeff_for_highT
         real(dp) :: logT_max_for_standard_mesh_delta_coeff
         real(dp) :: logT_min_for_highT_mesh_delta_coeff
         real(dp) :: remesh_dt_limit
         real(dp) :: mesh_Pgas_div_P_exponent

         logical :: okay_to_remesh, trace_mesh_adjust_error_in_conservation
         real(dp) :: max_rel_delta_IE_for_mesh_total_energy_balance
         integer :: max_allowed_nz

         logical :: restore_mesh_on_retry
         integer :: num_steps_to_hold_mesh_after_retry

         real(dp) :: mesh_max_allowed_ratio, remesh_max_allowed_logT

         real(dp) :: max_delta_x_for_merge

         logical :: show_mesh_changes

         logical :: mesh_adjust_use_quadratic, &
            mesh_adjust_get_T_from_E, mesh_ok_to_merge

         integer :: mesh_dump_call_number, mesh_max_k_old_for_split, mesh_min_k_old_for_split

         real(dp) :: max_dq, min_dq, min_dq_for_split, min_dq_for_logT, &
            min_dq_for_xa, min_dq_for_xa_convective

         real(dp) :: mesh_min_dlnR
         logical :: merge_if_dlnR_too_small

         real(dp) :: mesh_min_dr_div_dRstar
         logical :: merge_if_dr_div_dRstar_too_small

         real(dp) :: mesh_min_dr_div_cs
         logical :: merge_if_dr_div_cs_too_small

         real(dp) :: max_center_cell_dq
         real(dp) :: max_surface_cell_dq
         integer :: max_num_subcells
         integer :: max_num_merge_cells


         character (len=iso_name_length) :: mesh_logX_species(num_mesh_logX)
         real(dp) :: mesh_logX_min_for_extra(num_mesh_logX)
         real(dp) :: mesh_dlogX_dlogP_extra(num_mesh_logX)
         real(dp) :: mesh_dlogX_dlogP_full_on(num_mesh_logX)
         real(dp) :: mesh_dlogX_dlogP_full_off(num_mesh_logX)

         real(dp) :: mesh_dlog_eps_min_for_extra
         real(dp) :: mesh_dlog_eps_dlogP_full_on
         real(dp) :: mesh_dlog_eps_dlogP_full_off

         real(dp) :: mesh_dlog_pp_dlogP_extra
         real(dp) :: mesh_dlog_cno_dlogP_extra
         real(dp) :: mesh_dlog_3alf_dlogP_extra

         real(dp) :: mesh_dlog_burn_c_dlogP_extra
         real(dp) :: mesh_dlog_burn_n_dlogP_extra
         real(dp) :: mesh_dlog_burn_o_dlogP_extra
         real(dp) :: mesh_dlog_burn_ne_dlogP_extra
         real(dp) :: mesh_dlog_burn_na_dlogP_extra
         real(dp) :: mesh_dlog_burn_mg_dlogP_extra
         real(dp) :: mesh_dlog_burn_si_dlogP_extra
         real(dp) :: mesh_dlog_burn_s_dlogP_extra
         real(dp) :: mesh_dlog_burn_ar_dlogP_extra
         real(dp) :: mesh_dlog_burn_ca_dlogP_extra
         real(dp) :: mesh_dlog_burn_ti_dlogP_extra
         real(dp) :: mesh_dlog_burn_cr_dlogP_extra
         real(dp) :: mesh_dlog_burn_fe_dlogP_extra

         real(dp) :: mesh_dlog_cc_dlogP_extra
         real(dp) :: mesh_dlog_co_dlogP_extra
         real(dp) :: mesh_dlog_oo_dlogP_extra

         real(dp) :: mesh_dlog_pnhe4_dlogP_extra
         real(dp) :: mesh_dlog_photo_dlogP_extra
         real(dp) :: mesh_dlog_other_dlogP_extra
         
         integer :: mesh_delta_coeff_factor_smooth_iters

         real(dp) :: P_function_weight, E_function_weight, E_function_param

         real(dp) :: T_function1_weight
         real(dp) :: T_function2_weight
         real(dp) :: T_function2_param

         real(dp) :: R_function_weight
         real(dp) :: R_function_param

         real(dp) :: R_function2_weight
         real(dp) :: R_function2_param1
         real(dp) :: R_function2_param2

         real(dp) :: R_function3_weight

         real(dp) :: M_function_weight
         real(dp) :: M_function_param

         real(dp) :: gradT_function_weight
         real(dp) :: log_tau_function_weight
         real(dp) :: log_kap_function_weight
         real(dp) :: omega_function_weight

         real(dp) :: gam_function_weight, gam_function_param1, gam_function_param2

         character (len=iso_name_length) :: xa_function_species(num_xa_function)
         real(dp) :: xa_function_weight(num_xa_function)
         real(dp) :: xa_function_param(num_xa_function)
         real(dp) :: xa_mesh_delta_coeff(num_xa_function)
         
         logical :: use_split_merge_amr, trace_split_merge_amr, &
            equal_split_density_amr, split_merge_amr_skip_inner_cell, &
            split_merge_amr_log_zoning, &
            split_merge_amr_logtau_zoning, &
            split_merge_amr_okay_to_split_nz, &
            split_merge_amr_okay_to_split_1, &
            merge_amr_inhibit_at_jumps
         integer :: split_merge_amr_nz_baseline, split_merge_amr_max_iters
         real(dp) :: split_merge_amr_MaxLong, split_merge_amr_MaxShort, &
            merge_amr_max_abs_du_div_cs, split_merge_amr_dq_min, split_merge_amr_dq_max, &
            split_merge_amr_cntr_dq_min

         logical :: merge_amr_ignore_surface_cells
         integer :: merge_amr_k_for_ignore_surface_cells
         logical :: merge_amr_du_div_cs_limit_only_for_compression

         logical :: split_merge_amr_avoid_repeated_remesh

      ! nuclear reaction parameters

         character (len=strlen) :: screening_mode, default_net_name
         real(dp) :: net_logTcut_lo, net_logTcut_lim

         real(dp) :: reaction_neuQs_factor, nonlocal_NiCo_kap_gamma, &
            dtau_gamma_NiCo_decay_heat, max_logT_for_net
         
         logical :: nonlocal_NiCo_decay_heat

         real(dp) :: weak_rate_factor, mix_factor, dxdt_nuc_factor, &
            fe56ec_fake_factor, min_T_for_fe56ec_fake_factor, max_abar_for_burning


      ! extra heat near surface to model irradiation

         real(dp) :: column_depth_for_irradiation
         real(dp) :: irradiation_flux


      ! extra heat

         real(dp) :: inject_uniform_extra_heat, &
            min_q_for_uniform_extra_heat, &
            max_q_for_uniform_extra_heat

         real(dp) :: inject_extra_ergs_sec, &
            base_of_inject_extra_ergs_sec, &
            total_mass_for_inject_extra_ergs_sec, &
            start_time_for_inject_extra_ergs_sec, &
            duration_for_inject_extra_ergs_sec, &
            inject_until_reach_model_with_total_energy

      ! element diffusion parameters

         logical :: show_diffusion_info, show_diffusion_substep_info, &
            show_diffusion_timing

         logical :: do_element_diffusion
         logical :: diffusion_use_iben_macdonald
         logical :: diffusion_use_paquette
         logical :: diffusion_use_cgs_solver
         logical :: diffusion_use_full_net
         logical :: do_WD_sedimentation_heating
         logical :: do_diffusion_heating
         real(dp) :: min_xa_for_WD_sedimentation_heating
         real(dp) :: cgs_thermal_diffusion_eta_full_on
         real(dp) :: cgs_thermal_diffusion_eta_full_off
         real(dp) :: diffusion_dt_limit
         real(dp) :: diffusion_min_dq_at_surface, diffusion_dt_div_timescale, &
            diffusion_min_T_at_surface, diffusion_min_dq_ratio_at_surface
         integer :: diffusion_min_num_substeps

         real(dp) :: diffusion_min_X_hard_limit
         real(dp) :: diffusion_X_total_atol, diffusion_X_total_rtol
         real(dp) :: diffusion_v_max, diffusion_upwind_abs_v_limit
         real(dp) :: diffusion_gamma_full_off
         real(dp) :: diffusion_gamma_full_on
         real(dp) :: diffusion_T_full_on
         real(dp) :: diffusion_T_full_off
         real(dp) :: D_mix_ignore_diffusion

         real(dp) :: diffusion_Vlimit_dm_full_on
         real(dp) :: diffusion_Vlimit_dm_full_off
         real(dp) :: diffusion_Vlimit

         logical :: diffusion_calculates_ionization
         integer :: diffusion_nsmooth_typical_charge, &
            diffusion_max_iters_per_substep, &
            diffusion_max_retries_per_substep
         real(dp) :: &
            diffusion_tol_correction_max, diffusion_tol_correction_norm

         real(dp) :: diffusion_AD_dm_full_on
         real(dp) :: diffusion_AD_dm_full_off
         real(dp) :: diffusion_AD_boost_factor

         real(dp) :: diffusion_GT_factor, diffusion_SIG_factor

         real(dp) :: diffusion_max_T_for_radaccel, diffusion_min_T_for_radaccel
         integer :: diffusion_max_Z_for_radaccel, diffusion_min_Z_for_radaccel
         logical :: diffusion_screening_for_radaccel
         character (len=strlen) :: op_mono_data_path, op_mono_data_cache_filename

         integer :: diffusion_num_classes
         character (len=iso_name_length) :: &
            diffusion_class_representative(max_num_diffusion_classes)
         real(dp) :: diffusion_class_A_max(max_num_diffusion_classes)
         real(dp) :: diffusion_class_typical_charge(max_num_diffusion_classes)
         real(dp) :: diffusion_class_factor(max_num_diffusion_classes)

         logical :: diffusion_use_isolve
         real(dp) :: diffusion_rtol_for_isolve
         real(dp) :: diffusion_atol_for_isolve
         integer :: diffusion_maxsteps_for_isolve
         character(len=64) :: diffusion_isolve_solver


      ! eos controls
         logical :: &
            use_eosDT_ideal_gas, use_eosDT_HELMEOS, &
            report_eos_settings_at_start_of_run, &
            use_fixed_XZ_for_eos
         real(dp) :: fixed_X_for_eos, fixed_Z_for_eos

      ! opacity controls

         logical :: use_simple_es_for_kap, use_starting_composition_for_kap

         real(dp) :: min_kap_for_dPrad_dm_eqn
         real(dp) :: opacity_max, opacity_factor
         real(dp) :: min_logT_for_opacity_factor_off
         real(dp) :: min_logT_for_opacity_factor_on
         real(dp) :: max_logT_for_opacity_factor_on
         real(dp) :: max_logT_for_opacity_factor_off
         real(dp) :: &
            low_logT_op_mono_full_off, low_logT_op_mono_full_on, &
            high_logT_op_mono_full_off, high_logT_op_mono_full_on, &
            op_mono_min_X_to_include
         logical :: use_op_mono_alt_get_kap


      ! asteroseismology controls

         logical :: get_delta_nu_from_scaled_solar
         real(dp) :: nu_max_sun, delta_nu_sun, Teff_sun, delta_Pg_mode_freq


      ! hydro parameters

         real(dp) :: eps_nuc_factor, max_abs_eps_nuc, op_split_burn_eps_nuc_infall_limit
         real(dp) :: eps_WD_sedimentation_factor, non_nuc_neu_factor, eps_grav_factor, eps_mdot_factor
         real(dp) :: dedt_eqn_r_scale, &
            max_eta_for_dedt_form_of_energy_eqn, max_gamma_for_dedt_form_of_energy_eqn

         logical :: use_mass_corrections, use_gravity_rotation_correction, &
            use_dPrad_dm_form_of_T_gradient_eqn, &
            use_dedt_form_of_energy_eqn, &
            no_dedt_form_during_relax, &
            always_use_dedt_form_of_energy_eqn, &
            use_dedt_form_with_total_energy_conservation, &
            always_use_eps_grav_form_of_energy_eqn, &
            use_Fraley_time_centering, &
            include_P_in_Fraley_time_centering, &
            include_L_in_Fraley_time_centering, &
            use_time_centered_eps_grav
            
         integer :: steps_before_always_use_dedt_form_of_energy_eqn, &
            steps_before_use_Fraley_time_centering

         real(dp) :: &
            RTI_A, RTI_B, RTI_C, RTI_D, &
            RTI_C_X0_frac, RTI_C_X_factor, RTI_energy_floor, RTI_D_mix_floor, &
            RTI_max_alpha, RTI_min_dm_behind_shock_for_full_on, &
            RTI_dm_for_center_eta_nondecreasing, RTI_min_m_for_D_mix_floor, &
            RTI_log_max_boost, RTI_m_full_boost, RTI_m_no_boost

         logical :: include_composition_in_eps_grav
         real(dp) :: velocity_logT_lower_bound, &
            max_dt_yrs_for_velocity_logT_lower_bound, velocity_q_upper_bound
         real(dp) :: Gamma_lnS_eps_grav_full_on, Gamma_lnS_eps_grav_full_off
         integer :: max_num_surf_revisions
         real(dp) :: max_abs_rel_change_surf_lnS

         logical :: op_split_burn
         real(dp) :: op_split_burn_min_T, op_split_burn_eps, op_split_burn_odescal, &
            op_split_burn_min_T_for_variable_T_solver

         logical :: include_L_in_error_est, &
            include_v_in_error_est, include_u_in_error_est

         logical :: trace_solver_damping

         real(dp) :: tol_correction_norm
         real(dp) :: tol_max_correction
         real(dp) :: correction_xa_limit
         real(dp) :: tol_correction_high_T_limit, &
            tol_correction_norm_high_T, tol_max_correction_high_T
         real(dp) :: tol_correction_extreme_T_limit, &
            tol_correction_norm_extreme_T, tol_max_correction_extreme_T
         real(dp) :: tol_bad_max_correction
         integer :: bad_max_correction_series_limit

         real(dp) :: tol_residual_norm1, tol_max_residual1
         real(dp) :: tol_residual_norm2, tol_max_residual2
         real(dp) :: tol_residual_norm3, tol_max_residual3
         
         real(dp) :: warning_limit_for_max_residual

         logical :: relax_use_gold_tolerances
         real(dp) :: relax_tol_correction_norm, relax_tol_max_correction, &
            relax_tol_residual_norm1, relax_tol_max_residual1, &
            relax_tol_residual_norm2, relax_tol_max_residual2, &
            relax_tol_residual_norm3, relax_tol_max_residual3, &
            relax_maxT_for_gold_tolerances
         integer :: &
            relax_solver_iters_timestep_limit, &
            relax_iter_for_resid_tol2, relax_iter_for_resid_tol3
         
         integer :: steps_before_start_stress_test
         logical :: stress_test_relax
         
         logical :: use_gold_tolerances
         integer :: &
            gold_solver_iters_timestep_limit, &
            steps_before_use_gold_tolerances, &
            gold_iter_for_resid_tol2, &
            gold_iter_for_resid_tol3
         real(dp) :: &
            maxT_for_gold_tolerances, &
            gold_tol_residual_norm1, &
            gold_tol_max_residual1, &
            gold_tol_residual_norm2, &
            gold_tol_max_residual2, &
            gold_tol_residual_norm3, &
            gold_tol_max_residual3
         
         logical :: use_gold2_tolerances
         integer :: &
            gold2_solver_iters_timestep_limit, &
            steps_before_use_gold2_tolerances, &
            gold2_iter_for_resid_tol2, &
            gold2_iter_for_resid_tol3
         real(dp) :: &
            gold2_tol_residual_norm1, &
            gold2_tol_max_residual1, &
            gold2_tol_residual_norm2, &
            gold2_tol_max_residual2, &
            gold2_tol_residual_norm3, &
            gold2_tol_max_residual3

         logical :: convergence_ignore_equL_residuals, &
            convergence_ignore_alpha_RTI_residuals, &
            include_rotation_in_total_energy

         integer :: iter_for_resid_tol2, iter_for_resid_tol3

         integer :: solver_itermin, solver_itermin_until_reduce_min_corr_coeff
         logical :: do_solver_damping_for_neg_xa

         logical :: scale_max_correction_for_negative_surf_lum
         real(dp) :: max_frac_for_negative_surf_lum

         real(dp) :: &
            hydro_mtx_max_allowed_abs_dlogT, &
            hydro_mtx_max_allowed_abs_dlogRho, &
            hydro_mtx_max_allowed_abs_dlogPgas, &
            min_logT_for_hydro_mtx_max_allowed, &
            hydro_mtx_max_allowed_logT, &
            hydro_mtx_max_allowed_logRho, &
            hydro_mtx_max_allowed_logPgas, &
            hydro_mtx_min_allowed_logT, &
            hydro_mtx_min_allowed_logRho, &
            solver_reduced_min_corr_coeff, &
            solver_clip_dlogT, &
            solver_clip_dlogRho, &
            solver_clip_dlogPgas, &
            solver_clip_dlogR

         character (len=32) :: hydro_solver_name

         integer :: num_times_solver_reuse_mtx
         integer :: tiny_corr_coeff_limit
         real(dp) :: scale_correction_norm
         real(dp) :: corr_param_factor
         real(dp) :: scale_max_correction
         logical :: ignore_min_corr_coeff_for_scale_max_correction
         logical :: ignore_too_large_correction
         logical :: ignore_species_in_max_correction
         real(dp) :: corr_norm_jump_limit
         real(dp) :: max_corr_jump_limit
         real(dp) :: resid_norm_jump_limit
         real(dp) :: max_resid_jump_limit
         real(dp) :: corr_coeff_limit
         real(dp) :: tiny_corr_factor

         integer :: solver_max_tries_before_reject, max_tries1
         integer :: max_tries_for_retry, max_tries_after_5_retries, &
            max_tries_after_10_retries, max_tries_after_20_retries

         integer :: retry_limit, redo_limit

         logical :: use_avQ_art_visc
         real(dp) :: avQ_cq, avQ_zsh

         real(dp) :: sig_term_limit, sig_min_factor_for_high_Tcenter, &
            Tcenter_min_for_sig_min_factor_full_on, &
            Tcenter_max_for_sig_min_factor_full_off, &
            max_delta_m_to_bdy_for_sig_min_factor, &
            delta_m_lower_for_sig_min_factor, &
            delta_m_upper_for_sig_min_factor

         ! rotation
         logical :: set_uniform_am_nu_non_rot
         real(dp) :: uniform_am_nu_non_rot

         logical :: set_min_am_nu_non_rot
         real(dp) :: min_am_nu_non_rot, min_center_Ye_for_min_am_nu_non_rot

         logical :: set_min_D_mix, set_min_D_mix_below_Tmax, set_min_D_mix_in_H_He
         real(dp) :: min_D_mix, min_D_mix_below_Tmax, min_center_Ye_for_min_D_mix, &
            min_D_mix_in_H_He, mass_lower_limit_for_min_D_mix, &
            mass_upper_limit_for_min_D_mix

         real(dp) :: smooth_outer_xa_big, smooth_outer_xa_small

         real(dp) :: am_sig_term_limit
         real(dp) :: am_D_mix_factor, am_gradmu_factor, &
            am_nu_factor, am_nu_non_rotation_factor
         logical :: skip_rotation_in_convection_zones

         real(dp) :: am_nu_DSI_factor, am_nu_SH_factor, am_nu_SSI_factor, &
            am_nu_ES_factor, am_nu_GSF_factor, am_nu_ST_factor, am_nu_visc_factor

         real(dp) :: &
            am_nu_omega_rot_factor, am_nu_omega_non_rot_factor, &
            am_nu_j_rot_factor, am_nu_j_non_rot_factor

         real(dp) :: D_visc_factor, D_DSI_factor, D_SH_factor, D_SSI_factor, &
            D_ES_factor, D_GSF_factor, D_ST_factor

         integer :: smooth_D_DSI, smooth_D_SH, smooth_D_SSI, &
            smooth_D_ES, smooth_D_GSF, smooth_D_ST, smooth_nu_ST, smooth_D_omega, smooth_am_nu_rot

         real(dp) :: adjust_J_fraction, min_q_for_adjust_J_lost, &
            min_J_div_delta_J
         logical :: simple_i_rot_flag, do_adjust_J_lost, &
            premix_omega, recalc_mixing_info_each_substep
         real(dp) :: angular_momentum_error_warn, angular_momentum_error_retry
         integer :: max_mdot_redo_cnt
         real(dp) :: min_years_dt_for_redo_mdot, mdot_crit_boost_factor, &
            mdot_revise_factor, implicit_mdot_boost
         real(dp) :: surf_w_div_w_crit_limit, surf_w_div_w_crit_tol
         logical :: fitted_fp_ft_i_rot
         real(dp) :: w_div_wcrit_max, w_div_wcrit_max2
         real(dp) :: fp_min, ft_min, fp_error_limit, ft_error_limit

         real(dp) :: D_mix_rotation_max_logT_full_on, &
            D_mix_rotation_min_logT_full_off
         
         real(dp) :: D_omega_mixing_rate, max_q_for_D_omega_zero_in_convection_region
         logical :: D_omega_mixing_across_convection_boundary
         
         real(dp) :: nu_omega_mixing_rate, max_q_for_nu_omega_zero_in_convection_region
         logical :: nu_omega_mixing_across_convection_boundary

         real(dp) :: jacobian_clip_limit
         real(dp) :: sum_xa_tolerance, xa_clip_limit

         real(dp) :: min_xa_hard_limit, min_xa_hard_limit_for_highT
         real(dp) :: logT_max_for_min_xa_hard_limit, &
            logT_min_for_min_xa_hard_limit_for_highT

         real(dp) :: sum_xa_hard_limit, sum_xa_hard_limit_for_highT
         real(dp) :: logT_max_for_sum_xa_hard_limit, &
            logT_min_for_sum_xa_hard_limit_for_highT

         logical :: report_solver_progress
         logical :: report_ierr, report_bad_negative_xa

         logical :: solver_numerical_jacobian
         integer :: solver_jacobian_nzlo, solver_jacobian_nzhi
         logical :: solver_check_everything
         logical :: solver_inspect_soln_flag
         real(dp) :: solver_test_partials_dx_0
         integer :: solver_test_partials_k, solver_test_partials_call_number, &
            solver_test_partials_iter_number, solver_save_photo_call_number, &
            solver_test_partials_k_low, solver_test_partials_k_high
         character (len=strlen) :: solver_test_partials_var_name, solver_test_partials_equ_name, &
            solver_test_partials_sink_name, solver_test_partials_show_dx_var_name
         logical :: solver_test_partials_write_eos_call_info
         logical :: solver_test_eos_partials, solver_test_kap_partials, &
            solver_test_net_partials, solver_test_atm_partials

         logical :: solver_show_correction_info
         real(dp) :: solver_epsder_chem, solver_epsder_struct
         
         integer :: trace_k
         logical :: fill_arrays_with_NaNs
         logical :: zero_when_allocate, warn_when_get_a_bad_eos_result, &
            warn_rates_for_high_temp
         real(dp) :: warn_when_large_rel_run_E_err, &
            warn_when_large_virial_thm_rel_err, max_safe_logT_for_rates, &
            eps_mdot_leak_frac_factor
            
         logical :: stop_for_bad_nums
         integer :: diffusion_dump_call_number, energy_conservation_dump_model_number
         
      ! timestep parameters
      
         real(dp) :: time_delta_coeff ! like mesh_delta_coeff

         real(dp) :: max_timestep
         real(dp) :: max_years_for_timestep
         real(dp) :: hi_T_max_years_for_timestep
         real(dp) :: max_timestep_hi_T_limit

         real(dp) :: min_timestep_factor
         real(dp) :: max_timestep_factor
         real(dp) :: max_timestep_factor_at_high_T, &
            min_logT_for_max_timestep_factor_at_high_T
         real(dp) :: timestep_factor_for_retries
         integer :: retry_hold
         integer :: neg_mass_fraction_hold

         real(dp) :: timestep_dt_factor
         logical :: use_dt_low_pass_controller
         
         real(dp) :: force_timestep_min, &
            force_timestep_min_years, force_timestep_min_factor

         real(dp) :: varcontrol_target, varcontrol_dt_limit_ratio_hard_max, &
            min_allowed_varcontrol_target

         real(dp) :: xa_scale, abs_xa_corr_limit

         integer :: solver_iters_timestep_limit

         integer :: burn_steps_limit
         integer :: burn_steps_hard_limit

         integer :: diffusion_steps_limit
         integer :: diffusion_steps_hard_limit

         integer :: diffusion_iters_limit
         integer :: diffusion_iters_hard_limit

         real(dp) :: &
            dt_div_dt_cell_collapse_limit, &
            dt_div_dt_cell_collapse_hard_limit, &
            dt_div_min_dr_div_cs_limit, &
            dt_div_min_dr_div_cs_hard_limit, &
            min_abs_du_div_cs_for_dt_div_min_dr_div_cs_limit, &
            min_q_for_dt_div_min_dr_div_cs_limit, &
            max_q_for_dt_div_min_dr_div_cs_limit
         integer :: min_k_for_dt_div_min_dr_div_cs_limit

         real(dp) :: dX_mix_dist_limit

         real(dp) :: dH_limit_min_H
         real(dp) :: dH_limit
         real(dp) :: dH_hard_limit
         logical :: dH_decreases_only

         real(dp) :: dH_div_H_limit_min_H
         real(dp) :: dH_div_H_limit
         real(dp) :: dH_div_H_hard_limit

         real(dp) :: dHe_limit_min_He
         real(dp) :: dHe_limit
         real(dp) :: dHe_hard_limit
         logical :: dHe_decreases_only

         real(dp) :: dHe_div_He_limit_min_He
         real(dp) :: dHe_div_He_limit
         real(dp) :: dHe_div_He_hard_limit

         real(dp) :: dHe3_limit_min_He3
         real(dp) :: dHe3_limit
         real(dp) :: dHe3_hard_limit
         logical :: dHe3_decreases_only

         real(dp) :: dHe3_div_He3_limit_min_He3
         real(dp) :: dHe3_div_He3_limit
         real(dp) :: dHe3_div_He3_hard_limit

         real(dp) :: dX_limit_min_X
         real(dp) :: dX_limit
         real(dp) :: dX_hard_limit
         logical :: dX_decreases_only

         real(dp) :: dX_div_X_limit_min_X
         real(dp) :: dX_div_X_limit
         real(dp) :: dX_div_X_hard_limit
         real(dp) :: dX_div_X_at_high_T_limit
         real(dp) :: dX_div_X_at_high_T_hard_limit
         real(dp) :: dX_div_X_at_high_T_limit_lgT_min

         real(dp) :: dX_nuc_drop_min_X_limit
         real(dp) :: dX_nuc_drop_max_A_limit
         real(dp) :: dX_nuc_drop_limit, dX_nuc_drop_limit_at_high_T
         real(dp) :: dX_nuc_drop_hard_limit
         real(dp) :: dX_nuc_drop_min_yrs_for_dt

         real(dp) :: dL_div_L_limit_min_L
         real(dp) :: dL_div_L_limit
         real(dp) :: dL_div_L_hard_limit

         real(dp) :: delta_lgP_limit
         real(dp) :: delta_lgP_hard_limit, delta_lgP_limit_min_lgP

         real(dp) :: delta_lgRho_limit
         real(dp) :: delta_lgRho_hard_limit, delta_lgRho_limit_min_lgRho

         real(dp) :: delta_lgT_limit
         real(dp) :: delta_lgT_hard_limit, delta_lgT_limit_min_lgT

         real(dp) :: delta_lgE_limit
         real(dp) :: delta_lgE_hard_limit, delta_lgE_limit_min_lgE

         real(dp) :: delta_lgR_limit
         real(dp) :: delta_lgR_hard_limit, delta_lgR_limit_min_lgR

         real(dp) :: delta_Ye_highT_limit
         real(dp) :: delta_Ye_highT_hard_limit
         real(dp) :: minT_for_highT_Ye_limit

         real(dp) :: delta_lgL_nuc_limit
         real(dp) :: delta_lgL_nuc_hard_limit
         real(dp) :: delta_lgL_nuc_at_high_T_limit
         real(dp) :: delta_lgL_nuc_at_high_T_hard_limit
         real(dp) :: delta_lgL_nuc_at_high_T_limit_lgT_min
         real(dp) :: max_lgT_for_lgL_nuc_limit
         real(dp) :: lgL_nuc_burn_min, lgL_nuc_drop_factor
         
         real(dp) :: delta_lgL_power_photo_limit = -1
         real(dp) :: delta_lgL_power_photo_hard_limit = -1
         real(dp) :: min_lgT_for_lgL_power_photo_limit = 9d0
         real(dp) :: lgL_power_photo_burn_min = 10d0
         real(dp) :: lgL_power_photo_drop_factor = 10

         real(dp) :: delta_lgL_nuc_cat_limit
         real(dp) :: delta_lgL_nuc_cat_hard_limit
         real(dp) :: lgL_nuc_cat_burn_min
         real(dp) :: lgL_nuc_mix_dist_limit

         real(dp) :: delta_lgL_H_limit
         real(dp) :: delta_lgL_H_hard_limit, lgL_H_burn_relative_limit
         real(dp) :: lgL_H_burn_min, lgL_H_drop_factor

         real(dp) :: delta_lgL_He_limit
         real(dp) :: delta_lgL_He_hard_limit, lgL_He_burn_relative_limit
         real(dp) :: lgL_He_burn_min, lgL_He_drop_factor

         real(dp) :: delta_lgL_z_limit
         real(dp) :: delta_lgL_z_hard_limit, lgL_z_burn_relative_limit
         real(dp) :: lgL_z_burn_min, lgL_z_drop_factor

         real(dp) :: delta_lgP_cntr_limit
         real(dp) :: delta_lgP_cntr_hard_limit

         real(dp) :: delta_lgRho_cntr_limit
         real(dp) :: delta_lgRho_cntr_hard_limit

         real(dp) :: delta_lgT_cntr_limit
         real(dp) :: delta_lgT_cntr_hard_limit
         logical :: delta_lgT_cntr_limit_only_after_near_zams

         real(dp) :: delta_lgT_max_limit
         real(dp) :: delta_lgT_max_hard_limit
         real(dp) :: delta_lgT_max_limit_lgT_min
         logical :: delta_lgT_max_limit_only_after_near_zams

         real(dp) :: delta_lgT_max_at_high_T_limit
         real(dp) :: delta_lgT_max_at_high_T_hard_limit
         real(dp) :: delta_lgT_max_at_high_T_limit_lgT_min

         real(dp) :: delta_log_eps_nuc_limit
         real(dp) :: delta_log_eps_nuc_hard_limit
         
         logical :: delta_dX_div_X_drop_only
         logical :: delta_lg_XH_drop_only
         logical :: delta_lg_XHe_drop_only
         logical :: delta_lg_XC_drop_only
         logical :: delta_lg_XNe_drop_only
         logical :: delta_lg_XO_drop_only
         logical :: delta_lg_XSi_drop_only
         logical :: delta_XH_drop_only
         logical :: delta_XHe_drop_only
         logical :: delta_XC_drop_only
         logical :: delta_XNe_drop_only
         logical :: delta_XO_drop_only
         logical :: delta_XSi_drop_only

         real(dp) :: delta_dX_div_X_cntr_limit
         real(dp) :: delta_dX_div_X_cntr_hard_limit
         real(dp) :: delta_dX_div_X_cntr_min
         real(dp) :: delta_dX_div_X_cntr_max

         real(dp) :: delta_lg_XH_cntr_limit
         real(dp) :: delta_lg_XH_cntr_hard_limit
         real(dp) :: delta_lg_XH_cntr_min
         real(dp) :: delta_lg_XH_cntr_max

         real(dp) :: delta_lg_XHe_cntr_limit
         real(dp) :: delta_lg_XHe_cntr_hard_limit
         real(dp) :: delta_lg_XHe_cntr_min
         real(dp) :: delta_lg_XHe_cntr_max

         real(dp) :: delta_lg_XC_cntr_limit
         real(dp) :: delta_lg_XC_cntr_hard_limit
         real(dp) :: delta_lg_XC_cntr_min
         real(dp) :: delta_lg_XC_cntr_max

         real(dp) :: delta_lg_XNe_cntr_limit
         real(dp) :: delta_lg_XNe_cntr_hard_limit
         real(dp) :: delta_lg_XNe_cntr_min
         real(dp) :: delta_lg_XNe_cntr_max

         real(dp) :: delta_lg_XO_cntr_limit
         real(dp) :: delta_lg_XO_cntr_hard_limit
         real(dp) :: delta_lg_XO_cntr_min
         real(dp) :: delta_lg_XO_cntr_max

         real(dp) :: delta_lg_XSi_cntr_limit
         real(dp) :: delta_lg_XSi_cntr_hard_limit
         real(dp) :: delta_lg_XSi_cntr_min
         real(dp) :: delta_lg_XSi_cntr_max

         real(dp) :: delta_XH_cntr_limit
         real(dp) :: delta_XH_cntr_hard_limit
         real(dp) :: delta_XHe_cntr_limit
         real(dp) :: delta_XHe_cntr_hard_limit
         real(dp) :: delta_XC_cntr_limit
         real(dp) :: delta_XC_cntr_hard_limit
         real(dp) :: delta_XNe_cntr_limit
         real(dp) :: delta_XNe_cntr_hard_limit
         real(dp) :: delta_XO_cntr_limit
         real(dp) :: delta_XO_cntr_hard_limit
         real(dp) :: delta_XSi_cntr_limit
         real(dp) :: delta_XSi_cntr_hard_limit

         real(dp) :: delta_lgTeff_limit
         real(dp) :: delta_lgTeff_hard_limit

         real(dp) :: delta_lgL_limit, delta_lgL_limit_L_min
         real(dp) :: delta_lgL_hard_limit

         real(dp) :: delta_HR_ds_L
         real(dp) :: delta_HR_ds_Teff
         real(dp) :: delta_HR_limit
         real(dp) :: delta_HR_hard_limit

         real(dp) :: delta_mdot_atol, delta_mdot_rtol
         real(dp) :: delta_mdot_limit, delta_mdot_hard_limit

         real(dp) :: adjust_J_q_limit, adjust_J_q_hard_limit

         real(dp) :: delta_lg_star_mass_limit
         real(dp) :: delta_lg_star_mass_hard_limit
         real(dp) :: &
            limit_for_rel_error_in_energy_conservation, &
            hard_limit_for_rel_error_in_energy_conservation

         logical :: relax_hard_limits_after_retry, never_skip_hard_limits
         logical :: report_dt_hard_limit_retries, report_solver_dt_info
         logical :: trace_evolve, report_min_dr_div_cs

      ! misc

         real(dp) :: min_chem_eqn_scale

         character (len=strlen) :: zams_filename
         logical :: set_rho_to_dm_div_dV
         
         logical :: use_d_eos_dxa

         logical :: use_other_eos
         logical :: use_other_surface_PT
         logical :: use_other_kap
         logical :: use_other_diffusion
         logical :: use_other_diffusion_factor
         logical :: use_other_mlt
         logical :: use_other_adjust_mdot
         logical :: use_other_j_for_adjust_J_lost
         logical :: use_other_alpha_mlt
         logical :: use_other_am_mixing
         logical :: use_other_brunt, use_other_brunt_smoothing
         logical :: use_other_solver_monitor
         logical :: use_other_build_initial_model
         logical :: use_other_cgrav
         logical :: use_other_mesh_delta_coeff_factor
         logical :: use_other_momentum_implicit
         logical :: use_other_momentum
         logical :: use_other_energy_implicit
         logical :: use_other_energy
         logical :: use_other_pressure
         logical :: use_other_mesh_functions
         logical :: use_other_eps_grav
         logical :: use_other_gradr_factor
         logical :: use_other_D_mix
         logical :: use_other_neu
         logical :: use_other_net_get
         logical :: use_other_opacity_factor
         logical :: use_other_diffusion_coefficients
         logical :: use_other_pgstar_plots
         logical :: use_other_eval_fp_ft
         logical :: use_other_eval_i_rot
         logical :: use_other_remove_surface
         logical :: use_other_torque
         logical :: use_other_torque_implicit
         logical :: use_other_wind
         logical :: use_other_accreting_state
         logical :: use_other_after_struct_burn_mix
         logical :: use_other_before_struct_burn_mix
         logical :: use_other_export_pulse_data
         logical :: use_other_get_pulse_data
         logical :: use_other_edit_pulse_data
         logical :: use_other_astero_freq_corr
         logical :: use_other_timestep_limit

         real(dp) :: x_ctrl(num_x_ctrls)
         integer :: x_integer_ctrl(num_x_ctrls)
         logical :: x_logical_ctrl(num_x_ctrls)
         character (len=strlen) :: x_character_ctrl(num_x_ctrls)
