      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ! star_job controls

      character (len=1000) :: mesa_dir, &
         eosDT_cache_dir, &
         ionization_cache_dir, &
         kap_cache_dir, &
         rates_cache_dir

      logical :: pause_before_terminate

      character (len=1000) :: profile_columns_file, history_columns_file

      logical :: show_log_description_at_start, list_net_reactions
      logical :: show_net_reactions_info, show_net_species_info

      integer :: first_model_for_timing, steps_before_start_timing

      logical :: show_eqns_and_vars_names

      logical :: pgstar_flag, save_pgstar_files_when_terminate, &
         clear_initial_pgstar_history, clear_pgstar_history

      logical :: save_photo_when_terminate
      logical :: write_profile_when_terminate
      character (len=1000) :: filename_for_profile_when_terminate

      logical :: show_retry_counts_when_terminate, &
         show_timestep_limit_counts_when_terminate

      logical :: create_pre_main_sequence_model, pre_ms_relax_to_start_radiative_core
      real(dp) :: pre_ms_check_radiative_core_Lnuc_div_L_limit, pre_ms_check_radiative_core_min_mass
      integer :: pre_ms_relax_num_steps, pre_ms_min_steps_before_check_radiative_core
      real(dp) :: pre_ms_T_c, pre_ms_guess_rho_c, pre_ms_d_log10_P
      real(dp) :: pre_ms_logT_surf_limit, pre_ms_logP_surf_limit
      real(dp) :: pre_ms_check_radiative_core_start, pre_ms_check_radiative_core_stop

      logical :: create_initial_model
      integer :: initial_model_relax_num_steps, max_tries_for_create_initial_model
      real(dp) :: &
         radius_in_cm_for_create_initial_model, &
         mass_in_gm_for_create_initial_model, &
         center_logP_1st_try_for_create_initial_model, &
         entropy_1st_try_for_create_initial_model, &
         abs_e01_tolerance_for_create_initial_model, &
         abs_e02_tolerance_for_create_initial_model
      real(dp) :: initial_model_eps

      logical :: save_star_job_namelist
      character (len=strlen) :: star_job_namelist_name
      
      character (len=strlen) :: echo_at_start, echo_at_end

      logical :: load_saved_model, load_saved_model_for_RSP
      character (len=strlen) :: saved_model_name

      logical :: create_merger_model
      character (len=strlen) :: saved_model_for_merger_1
      character (len=strlen) :: saved_model_for_merger_2

      logical :: set_max_dt_to_frac_lifetime
      real(dp) :: max_frac_of_lifetime_per_step

      logical :: astero_just_call_my_extras_check_model

      logical :: relax_mass, relax_initial_mass
      real(dp) :: new_mass
      real(dp) :: lg_max_abs_mdot
      logical :: relax_mass_scale, relax_initial_mass_scale
      real(dp) :: dlgm_per_step, relax_M_center_dt, change_mass_years_for_dt
      logical :: relax_M_center, relax_initial_M_center
     
      logical :: relax_mass_to_remove_H_env, relax_initial_mass_to_remove_H_env
      real(dp) :: extra_mass_retained_by_remove_H_env

      logical :: relax_core, relax_initial_core
      real(dp) :: new_core_mass, dlg_core_mass_per_step, &
         relax_core_years_for_dt, core_avg_rho, core_avg_eps

      logical :: relax_R_center, relax_initial_R_center
      real(dp) :: new_R_center
      real(dp) :: dlgR_per_step, relax_R_center_dt

      logical :: zero_alpha_RTI, zero_initial_alpha_RTI

      logical :: set_v_center, set_initial_v_center
      logical :: relax_v_center, relax_initial_v_center
      real(dp) :: new_v_center
      real(dp) :: dv_per_step, relax_v_center_dt

      logical :: set_L_center, set_initial_L_center
      logical :: relax_L_center, relax_initial_L_center
      real(dp) :: new_L_center
      real(dp) :: dlgL_per_step, relax_L_center_dt

      integer :: remove_center_at_cell_k, remove_initial_center_at_cell_k
      logical :: remove_center_at_max_infall_speed, remove_initial_center_at_cmax_infall_speed
      real(dp) :: remove_center_by_temperature, remove_initial_center_by_temperature
      real(dp) :: remove_center_by_mass_fraction_q, remove_initial_center_by_mass_fraction_q
      real(dp) :: remove_center_by_delta_mass_gm, remove_initial_center_by_delta_mass_gm
      real(dp) :: remove_center_by_delta_mass_Msun, remove_initial_center_by_delta_mass_Msun
      real(dp) :: remove_center_by_mass_gm, remove_initial_center_by_mass_gm
      real(dp) :: remove_center_by_mass_Msun, remove_initial_center_by_mass_Msun
      real(dp) :: remove_center_by_radius_cm, remove_initial_center_by_radius_cm
      real(dp) :: remove_center_by_radius_Rsun, remove_initial_center_by_radius_Rsun
      real(dp) :: remove_center_by_he4, remove_initial_center_by_he4
      real(dp) :: remove_center_by_si28, remove_initial_center_by_si28
      real(dp) :: remove_center_to_reduce_co56_ni56, remove_initial_center_to_reduce_co56_ni56
      real(dp) :: remove_center_by_ye, remove_initial_center_by_ye
      real(dp) :: remove_center_by_entropy, remove_initial_center_by_entropy
      real(dp) :: remove_center_by_infall_kms, remove_initial_center_by_infall_kms
      real(dp) :: remove_center_logRho_limit, remove_fallback_speed_limit
      logical :: remove_initial_center_at_inner_max_abs_v, remove_center_at_inner_max_abs_v, &
         remove_fallback_at_each_step, remove_center_adjust_L_center, &
         remove_center_set_zero_v_center, remove_initial_fe_core, remove_fe_core, &
         retain_fallback_at_each_step, fallback_check_total_energy
         
      real(dp) :: limit_center_logP_at_each_step
      real(dp) :: min_q_for_remove_surface_by_v_surf_div_v_escape, &
         max_q_for_remove_surface_by_v_surf_div_v_escape

      integer :: remove_surface_at_cell_k, remove_initial_surface_at_cell_k
      real(dp) :: remove_surface_at_he_core_boundary, remove_initial_surface_at_he_core_boundary
      real(dp) :: remove_surface_by_optical_depth, remove_initial_surface_by_optical_depth
      real(dp) :: remove_surface_by_density, remove_initial_surface_by_density
      real(dp) :: remove_surface_by_pressure, remove_initial_surface_by_pressure
      real(dp) :: remove_surface_by_mass_fraction_q, remove_initial_surface_by_mass_fraction_q
      real(dp) :: remove_surface_by_mass_gm, remove_initial_surface_by_mass_gm
      real(dp) :: remove_surface_by_radius_cm, remove_initial_surface_by_radius_cm
      real(dp) :: remove_surface_by_mass_Msun, remove_initial_surface_by_mass_Msun
      real(dp) :: remove_surface_by_radius_Rsun, remove_initial_surface_by_radius_Rsun
      real(dp) :: remove_surface_by_v_surf_km_s, remove_initial_surface_by_v_surf_km_s
      real(dp) :: remove_surface_by_v_surf_div_v_escape, remove_initial_surface_by_v_surf_div_v_escape
      real(dp) :: remove_surface_by_v_surf_div_cs, remove_initial_surface_by_v_surf_div_cs

      real(dp) :: zero_initial_inner_v_by_mass_Msun, &
         zero_inner_v_by_mass_Msun

      logical :: report_mass_not_fe56

      logical :: relax_dxdt_nuc_factor, relax_initial_dxdt_nuc_factor
      real(dp) :: new_dxdt_nuc_factor
      real(dp) :: dxdt_nuc_factor_multiplier

      logical :: relax_eps_nuc_factor, relax_initial_eps_nuc_factor
      real(dp) :: new_eps_nuc_factor
      real(dp) :: eps_nuc_factor_multiplier

      logical :: relax_opacity_max, relax_initial_opacity_max
      real(dp) :: new_opacity_max
      real(dp) :: opacity_max_multiplier

      logical :: relax_max_surf_dq, relax_initial_max_surf_dq
      real(dp) :: new_max_surf_dq
      real(dp) :: max_surf_dq_multiplier

      logical :: relax_opacity_factor, relax_initial_opacity_factor
      real(dp) :: relax_to_this_opacity_factor, d_opacity_factor

      logical :: relax_tau_factor, relax_initial_tau_factor
      logical :: set_tau_factor, set_initial_tau_factor
      real(dp) :: relax_to_this_tau_factor, set_to_this_tau_factor
      real(dp) :: dlogtau_factor, &
         set_tau_factor_after_core_He_burn, set_tau_factor_after_core_C_burn, &
         relax_tau_factor_after_core_He_burn, relax_tau_factor_after_core_C_burn

      logical :: adjust_tau_factor_to_surf_density
      real(dp) :: base_for_adjust_tau_factor_to_surf_density

      logical :: relax_Tsurf_factor, relax_initial_Tsurf_factor
      logical :: set_Tsurf_factor, set_initial_Tsurf_factor
      real(dp) :: relax_to_this_Tsurf_factor, set_to_this_Tsurf_factor
      real(dp) :: dlogTsurf_factor

      logical :: relax_irradiation, relax_initial_irradiation
      logical :: set_irradiation, set_initial_irradiation
      integer :: relax_irradiation_min_steps
      real(dp) :: relax_to_this_irrad_flux, set_to_this_irrad_flux, &
         irrad_col_depth, relax_irradiation_max_yrs_dt

      logical :: relax_mass_change, relax_initial_mass_change
      integer :: relax_mass_change_min_steps
      real(dp) :: relax_mass_change_max_yrs_dt, &
         relax_mass_change_init_mdot, relax_mass_change_final_mdot

      logical :: change_RTI_flag, change_initial_RTI_flag, new_RTI_flag
      logical :: create_RSP_model
      logical :: change_RSP_flag, change_initial_RSP_flag, new_RSP_flag
      logical :: change_v_flag, change_initial_v_flag, new_v_flag
      logical :: change_D_omega_flag, change_initial_D_omega_flag, new_D_omega_flag
      logical :: change_am_nu_rot_flag, change_initial_am_nu_rot_flag, &
         new_am_nu_rot_flag, use_D_omega_for_am_nu_rot
      real(dp) :: center_ye_limit_for_v_flag
      logical :: change_u_flag, change_initial_u_flag, new_u_flag
      logical :: change_reconstruction_flag, change_initial_reconstruction_flag, &
         new_reconstruction_flag

      logical :: change_rotation_flag, change_initial_rotation_flag, new_rotation_flag

      logical :: set_omega, set_initial_omega
      integer :: set_omega_step_limit, set_near_zams_omega_steps
      real(dp) :: new_omega

      logical :: set_omega_div_omega_crit, set_initial_omega_div_omega_crit
      integer :: set_omega_div_omega_crit_step_limit, set_near_zams_omega_div_omega_crit_steps
      real(dp) :: new_omega_div_omega_crit

      logical :: set_surface_rotation_v, set_initial_surface_rotation_v
      integer :: set_surf_rotation_v_step_limit, set_near_zams_surface_rotation_v_steps
      real(dp) :: new_surface_rotation_v

      logical :: &
         relax_omega, &
         relax_initial_omega, &
         near_zams_relax_omega, &
         relax_omega_div_omega_crit, &
         relax_initial_omega_div_omega_crit, &
         near_zams_relax_omega_div_omega_crit, &
         relax_surface_rotation_v, &
         relax_initial_surface_rotation_v, &
         near_zams_relax_initial_surface_rotation_v
      integer :: num_steps_to_relax_rotation
      real(dp) :: relax_omega_max_yrs_dt

      logical :: set_uniform_initial_composition
      real(dp) :: initial_h1, initial_h2, initial_he3, initial_he4
      integer :: initial_zfracs
      logical :: dump_missing_metals_into_heaviest

      logical :: set_abundance, set_initial_abundance
      character (len=5) :: chem_name
      real(dp) :: new_frac
      integer :: set_abundance_nzlo, set_abundance_nzhi

      logical :: replace_element, replace_initial_element
      character (len=5) :: chem_name1, chem_name2
      integer :: replace_element_nzlo, replace_element_nzhi

      logical :: relax_initial_composition, relax_initial_to_xaccrete
      character (len=strlen) :: relax_composition_filename
      integer :: num_steps_to_relax_composition
      real(dp) :: timescale_for_relax_composition

      real(dp) :: report_cell_for_xm
      logical :: set_to_xa_for_accretion, set_initial_to_xa_for_accretion
      integer :: set_nzlo, set_nzhi

      logical :: change_Y, change_initial_Y
      logical :: relax_Y, relax_initial_Y
      real(dp) :: relax_Y_minq, relax_Y_maxq, new_Y

      logical :: change_Z, change_initial_Z
      logical :: relax_Z, relax_initial_Z
      real(dp) :: relax_Z_minq, relax_Z_maxq, new_Z

      logical :: relax_initial_angular_momentum
      character (len=strlen) :: relax_angular_momentum_filename
      integer :: max_steps_to_relax_angular_momentum
      real(dp) :: timescale_for_relax_angular_momentum
      real(dp) :: max_dt_for_relax_angular_momentum
      real(dp) :: num_timescales_for_relax_angular_momentum
      logical :: relax_angular_momentum_constant_omega_center

      logical :: relax_initial_entropy
      integer :: max_steps_to_relax_entropy
      real(dp) :: timescale_for_relax_entropy
      real(dp) :: max_dt_for_relax_entropy
      real(dp) :: num_timescales_for_relax_entropy
      character (len=strlen) :: relax_entropy_filename
      character (len=strlen) :: get_entropy_for_relax_from_eos

      integer :: steps_to_take_before_terminate
      character (len=strlen) :: stop_if_this_file_exists

      logical :: set_initial_age
      real(dp) :: initial_age

      logical :: set_initial_model_number
      integer :: initial_model_number

      logical :: set_initial_number_retries
      integer :: initial_number_retries
      
      logical :: set_initial_dt, limit_initial_dt
      real(dp) :: years_for_initial_dt, seconds_for_initial_dt

      logical :: set_HELM_OPAL_Zs
      real(dp) :: Z_all_HELM
      real(dp) :: Z_all_OPAL

      logical :: set_initial_cumulative_energy_error, set_cumulative_energy_error, &
         set_cumulative_energy_error_each_relax
      integer :: set_cumulative_energy_error_at_step
      real(dp) :: new_cumulative_energy_error, &
         set_cumulative_energy_error_each_step_if_age_less_than

      logical :: change_net, change_initial_net
      logical :: change_small_net, change_initial_small_net
      logical :: adjust_abundances_for_new_isos
      character (len=strlen) :: new_net_name, new_small_net_name, &
         h_he_net, co_net, adv_net

      logical :: set_rates_preference
      integer :: new_rates_preference
      character (len=16) :: set_rate_c12ag, set_rate_n14pg, set_rate_3a, set_rate_1212

      logical :: set_uniform_xa_from_file, set_uniform_initial_xa_from_file
      character (len=strlen) :: file_for_uniform_xa

      logical :: mix_section, mix_initial_section
      integer :: mix_section_nzlo, mix_section_nzhi

      integer :: weaklib_blend_hi_Z
      real(dp) :: &
         T9_weaklib_full_off, T9_weaklib_full_on, &
         T9_weaklib_full_off_hi_Z, T9_weaklib_full_on_hi_Z

      logical :: use_suzuki_weak_rates

      logical :: use_special_weak_rates
      character (len=1000) :: special_weak_states_file, special_weak_transitions_file
      character (len=strlen) :: ion_coulomb_corrections, electron_coulomb_corrections

      real(dp) :: mix_envelope_down_to_T, mix_initial_envelope_down_to_T

      logical :: auto_extend_net

      logical :: enable_adaptive_network
      real(dp) :: min_x_for_keep, min_x_for_n, min_x_for_add, &
         max_Z_for_add, max_N_for_add, max_A_for_add

      integer :: save_model_number
      character (len=strlen) :: save_model_filename
      logical :: save_model_when_terminate

      character (len=strlen) :: required_termination_code_string(num_termination_code_strings)

      logical :: profile_starting_model
      integer :: profile_model_number

      integer :: internals_num

      logical :: report_retries

      character (len=strlen) :: net_reaction_filename, jina_reaclib_filename
      real(dp) :: jina_reaclib_min_T9

      character (len=strlen) :: rate_tables_dir, rate_cache_suffix

      logical :: read_extra_star_job_inlist1
      character (len=strlen) :: extra_star_job_inlist1_name

      logical :: read_extra_star_job_inlist2
      character (len=strlen) :: extra_star_job_inlist2_name

      logical :: read_extra_star_job_inlist3
      character (len=strlen) :: extra_star_job_inlist3_name

      logical :: read_extra_star_job_inlist4
      character (len=strlen) :: extra_star_job_inlist4_name

      logical :: read_extra_star_job_inlist5
      character (len=strlen) :: extra_star_job_inlist5_name

      integer :: save_pulse_data_for_model_number
      logical :: save_pulse_data_when_terminate
      character (len=strlen) :: save_pulse_data_filename

      character (len=strlen) :: chem_isotopes_filename

      character(len=strlen) :: eos_file_prefix
      character(len=strlen) :: ionization_file_prefix, ionization_Z1_suffix

      integer :: extras_lipar, extras_lrpar, extras_lcpar, extras_llpar
      integer :: extras_ipar(max_extras_params)
      real(dp) :: extras_rpar(max_extras_params)
      character(len=max_extras_cpar_len) :: extras_cpar(max_extras_params)
      logical :: extras_lpar(max_extras_params)

      integer :: num_special_rate_factors
      real(dp) :: special_rate_factor(max_num_special_rate_factors)
      character(len=maxlen_reaction_Name) :: &
         reaction_for_special_factor(max_num_special_rate_factors)

      logical :: do_special_test

      integer :: color_num_files
      character(len=strlen) :: color_file_names(max_num_color_files)
      integer :: color_num_colors(max_num_bcs_per_file)


      logical :: warn_run_star_extras
      
      ! When to run cleanup to deallocate unneeded data
      logical :: report_garbage_collection
      integer :: num_steps_for_garbage_collection
      
