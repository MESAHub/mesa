
      character(len=24) :: version_number ! mesa version from file $MESA_DIR/data/version_number
      
      integer, parameter :: star_def_version = 17
      
      integer, parameter :: nz_alloc_extra = 200
      
      
      integer, parameter :: XL = dp
      !integer, parameter :: XL = qp ! will need some code changes to use quad precision

      ! marsaglia and zaman random number generator. period is 2**43 with
      ! 900 million different sequences. the state of the generator (for restarts)
      integer, parameter :: rand_u_len=97

      ! result codes for various evolution routines
      integer, parameter :: keep_going = 0
      integer, parameter :: redo = 1 
         ! repeat current step with same timestep
         ! for example, after changing a parameter such as mdot
      integer, parameter :: retry = 2
         ! repeat current step with smaller timestep
      integer, parameter :: terminate = 3

      ! phases of evolution
      integer, parameter :: phase_relax = -1
      integer, parameter :: phase_starting = 1
      integer, parameter :: phase_PreMS = phase_starting + 1
      integer, parameter :: phase_ZAMS = phase_PreMS + 1
      integer, parameter :: phase_IAMS = phase_ZAMS + 1
      integer, parameter :: phase_TAMS = phase_IAMS + 1
      integer, parameter :: phase_He_Burn = phase_TAMS + 1
      integer, parameter :: phase_ZACHeB = phase_He_Burn + 1
      integer, parameter :: phase_TACHeB = phase_ZACHeB + 1
      integer, parameter :: phase_TP_AGB = phase_TACHeB + 1
      integer, parameter :: phase_C_Burn = phase_TP_AGB + 1
      integer, parameter :: phase_Ne_Burn = phase_C_Burn + 1
      integer, parameter :: phase_O_Burn = phase_Ne_Burn + 1
      integer, parameter :: phase_Si_Burn = phase_O_Burn + 1
      integer, parameter :: phase_WDCS = phase_Si_Burn + 1
      
      ! time_step limit identifiers
      integer, parameter :: Tlim_struc = 1
      integer, parameter :: Tlim_max_timestep_factor = Tlim_struc + 1
      integer, parameter :: Tlim_min_timestep_factor = Tlim_max_timestep_factor + 1
      integer, parameter :: Tlim_solver = Tlim_min_timestep_factor + 1
      integer, parameter :: Tlim_num_burn_steps = Tlim_solver + 1
      integer, parameter :: Tlim_num_diff_solver_steps = Tlim_num_burn_steps + 1
      integer, parameter :: Tlim_dX = Tlim_num_diff_solver_steps + 1
      integer, parameter :: Tlim_dH = Tlim_dX + 1
      integer, parameter :: Tlim_dHe = Tlim_dH + 1
      integer, parameter :: Tlim_dH_div_H = Tlim_dHe + 1
      integer, parameter :: Tlim_dHe_div_He = Tlim_dH_div_H + 1
      integer, parameter :: Tlim_dX_div_X = Tlim_dHe_div_He + 1
      integer, parameter :: Tlim_dL_div_L = Tlim_dX_div_X + 1
      integer, parameter :: Tlim_dlgP = Tlim_dL_div_L + 1
      integer, parameter :: Tlim_dlgRho = Tlim_dlgP + 1
      integer, parameter :: Tlim_dlgE = Tlim_dlgRho + 1
      integer, parameter :: Tlim_dlgT = Tlim_dlgE + 1
      integer, parameter :: Tlim_dlgR = Tlim_dlgT + 1
      integer, parameter :: Tlim_dlgL_nuc_cat = Tlim_dlgR + 1
      integer, parameter :: Tlim_dlgL_H = Tlim_dlgL_nuc_cat + 1
      integer, parameter :: Tlim_dlgL_He = Tlim_dlgL_H + 1
      integer, parameter :: Tlim_dlgL_z = Tlim_dlgL_He + 1
      integer, parameter :: Tlim_dlgL_nuc = Tlim_dlgL_z + 1
      integer, parameter :: Tlim_dlgTeff = Tlim_dlgL_nuc + 1
      integer, parameter :: Tlim_dlgRho_cntr = Tlim_dlgTeff + 1
      integer, parameter :: Tlim_dlgP_cntr = Tlim_dlgRho_cntr + 1
      integer, parameter :: Tlim_dlgT_cntr = Tlim_dlgP_cntr + 1
      integer, parameter :: Tlim_dlgT_max = Tlim_dlgT_cntr + 1
      integer, parameter :: Tlim_dlgT_max_at_high_T = Tlim_dlgT_max + 1
      integer, parameter :: Tlim_lg_XH_cntr = Tlim_dlgT_max_at_high_T + 1
      integer, parameter :: Tlim_dmstar = Tlim_lg_XH_cntr + 1
      integer, parameter :: Tlim_dt_div_dt_cell_collapse = Tlim_dmstar + 1
      integer, parameter :: Tlim_dt_div_min_dr_div_cs = Tlim_dt_div_dt_cell_collapse + 1
      integer, parameter :: Tlim_lgL = Tlim_dt_div_min_dr_div_cs + 1
      integer, parameter :: Tlim_force_timestep = Tlim_lgL + 1
      integer, parameter :: Tlim_max_timestep = Tlim_force_timestep + 1
      integer, parameter :: Tlim_timestep_hold = Tlim_max_timestep + 1
      integer, parameter :: Tlim_dX_div_X_cntr = Tlim_timestep_hold + 1
      integer, parameter :: Tlim_lg_XHe_cntr = Tlim_dX_div_X_cntr + 1
      integer, parameter :: Tlim_lg_XC_cntr = Tlim_lg_XHe_cntr + 1
      integer, parameter :: Tlim_lg_XNe_cntr = Tlim_lg_XC_cntr + 1
      integer, parameter :: Tlim_lg_XO_cntr = Tlim_lg_XNe_cntr + 1
      integer, parameter :: Tlim_lg_XSi_cntr = Tlim_lg_XC_cntr + 1
      integer, parameter :: Tlim_XH_cntr = Tlim_lg_XSi_cntr + 1
      integer, parameter :: Tlim_XHe_cntr = Tlim_XH_cntr + 1
      integer, parameter :: Tlim_XC_cntr = Tlim_XHe_cntr + 1
      integer, parameter :: Tlim_XNe_cntr = Tlim_XC_cntr + 1
      integer, parameter :: Tlim_XO_cntr = Tlim_XNe_cntr + 1
      integer, parameter :: Tlim_XSi_cntr = Tlim_XO_cntr + 1
      integer, parameter :: Tlim_dX_nuc_drop = Tlim_XSi_cntr + 1
      integer, parameter :: Tlim_neg_X = Tlim_dX_nuc_drop + 1
      integer, parameter :: Tlim_bad_Xsum = Tlim_neg_X + 1
      integer, parameter :: Tlim_dlgL_power_photo = Tlim_bad_Xsum + 1
      integer, parameter :: Tlim_delta_HR = Tlim_dlgL_power_photo + 1
      integer, parameter :: Tlim_dHe3 = Tlim_delta_HR + 1
      integer, parameter :: Tlim_dHe3_div_He3 = Tlim_dHe3 + 1
      integer, parameter :: Tlim_dlog_eps_nuc = Tlim_dHe3_div_He3 + 1
      integer, parameter :: Tlim_num_diff_solver_iters = Tlim_dlog_eps_nuc + 1
      integer, parameter :: Tlim_del_mdot = Tlim_num_diff_solver_iters + 1
      integer, parameter :: Tlim_adjust_J_q = Tlim_del_mdot + 1
      integer, parameter :: Tlim_delta_Ye_highT = Tlim_adjust_J_q + 1
      integer, parameter :: Tlim_error_in_energy_conservation = Tlim_delta_Ye_highT + 1
      integer, parameter :: Tlim_retry = Tlim_error_in_energy_conservation + 1
      integer, parameter :: Tlim_binary = Tlim_retry + 1
      integer, parameter :: Tlim_error_other = Tlim_binary + 1
      integer, parameter :: Tlim_other_timestep_limit = Tlim_error_other + 1
      
      
      integer, parameter :: numTlim = Tlim_other_timestep_limit
      
      character (len=24) :: dt_why_str(numTlim) ! indicates the reason for the timestep choice
      
   
      
      ! termination codes
      integer, parameter :: t_max_age = 1
      integer, parameter :: t_max_omega_div_omega_crit = t_max_age + 1
      integer, parameter :: t_peak_burn_vconv_div_cs_limit = t_max_omega_div_omega_crit + 1
      integer, parameter :: t_max_model_number = t_peak_burn_vconv_div_cs_limit + 1 
      integer, parameter :: t_eta_center_limit = t_max_model_number + 1 
      integer, parameter :: t_log_center_temp_limit = t_eta_center_limit + 1 
      integer, parameter :: t_log_center_temp_lower_limit = t_log_center_temp_limit + 1 
      integer, parameter :: t_log_center_density_lower_limit = t_log_center_temp_lower_limit + 1 
      integer, parameter :: t_log_center_density_limit = t_log_center_density_lower_limit + 1 
      integer, parameter :: t_gamma_center_limit = t_log_center_density_limit + 1 
      integer, parameter :: t_log_max_temp_upper_limit = t_gamma_center_limit + 1 
      integer, parameter :: t_log_max_temp_lower_limit = t_log_max_temp_upper_limit + 1 
      integer, parameter :: t_HB_limit = t_log_max_temp_lower_limit + 1 
      integer, parameter :: t_star_mass_min_limit = t_HB_limit + 1 
      integer, parameter :: t_star_mass_max_limit = t_star_mass_min_limit + 1 
      integer, parameter :: t_remnant_mass_min_limit = t_star_mass_max_limit + 1 
      integer, parameter :: t_ejecta_mass_max_limit = t_remnant_mass_min_limit + 1 
      integer, parameter :: t_star_species_mass_min_limit = t_ejecta_mass_max_limit + 1 
      integer, parameter :: t_star_species_mass_max_limit = t_star_species_mass_min_limit + 1       
      integer, parameter :: t_xmstar_min_limit = t_star_species_mass_max_limit + 1 
      integer, parameter :: t_xmstar_max_limit = t_xmstar_min_limit + 1 
      integer, parameter :: t_envelope_mass_limit = t_xmstar_max_limit + 1 
      integer, parameter :: t_envelope_fraction_left_limit = t_envelope_mass_limit + 1 
      integer, parameter :: t_he_core_mass_limit = t_envelope_fraction_left_limit + 1 
      integer, parameter :: t_co_core_mass_limit = t_he_core_mass_limit + 1 
      integer, parameter :: t_one_core_mass_limit = t_co_core_mass_limit + 1 
      integer, parameter :: t_fe_core_mass_limit = t_one_core_mass_limit + 1 
      integer, parameter :: t_neutron_rich_core_mass_limit = t_fe_core_mass_limit + 1 
      integer, parameter :: t_he_layer_mass_lower_limit = t_neutron_rich_core_mass_limit + 1 
      integer, parameter :: t_abs_diff_lg_LH_lg_Ls_limit = t_he_layer_mass_lower_limit + 1 
      integer, parameter :: t_Teff_lower_limit = t_abs_diff_lg_LH_lg_Ls_limit + 1 
      integer, parameter :: t_Teff_upper_limit = t_Teff_lower_limit + 1 
      integer, parameter :: t_delta_nu_lower_limit = t_Teff_upper_limit + 1 
      integer, parameter :: t_delta_nu_upper_limit = t_delta_nu_lower_limit + 1 
      integer, parameter :: t_delta_Pg_lower_limit = t_delta_nu_upper_limit + 1 
      integer, parameter :: t_delta_Pg_upper_limit = t_delta_Pg_lower_limit + 1 
      integer, parameter :: t_shock_mass_upper_limit = t_delta_Pg_upper_limit + 1 
      integer, parameter :: t_photosphere_m_sub_M_center_limit = t_shock_mass_upper_limit + 1 
      integer, parameter :: t_photosphere_m_lower_limit = t_photosphere_m_sub_M_center_limit + 1 
      integer, parameter :: t_photosphere_m_upper_limit = t_photosphere_m_lower_limit + 1 
      integer, parameter :: t_photosphere_r_lower_limit = t_photosphere_m_upper_limit + 1 
      integer, parameter :: t_photosphere_r_upper_limit = t_photosphere_r_lower_limit + 1 
      integer, parameter :: t_log_Teff_lower_limit = t_photosphere_r_upper_limit + 1 
      integer, parameter :: t_log_Teff_upper_limit = t_log_Teff_lower_limit + 1 
      integer, parameter :: t_log_Tsurf_lower_limit = t_log_Teff_upper_limit + 1 
      integer, parameter :: t_log_Tsurf_upper_limit = t_log_Tsurf_lower_limit + 1 
      integer, parameter :: t_log_Rsurf_lower_limit = t_log_Tsurf_upper_limit + 1 
      integer, parameter :: t_log_Rsurf_upper_limit = t_log_Rsurf_lower_limit + 1 
      integer, parameter :: t_log_Psurf_lower_limit = t_log_Rsurf_upper_limit + 1 
      integer, parameter :: t_log_Psurf_upper_limit = t_log_Psurf_lower_limit + 1 
      integer, parameter :: t_log_Dsurf_lower_limit = t_log_Psurf_upper_limit + 1 
      integer, parameter :: t_log_Dsurf_upper_limit = t_log_Dsurf_lower_limit + 1 
      integer, parameter :: t_log_L_lower_limit = t_log_Dsurf_upper_limit + 1 
      integer, parameter :: t_log_L_upper_limit = t_log_L_lower_limit + 1 
      integer, parameter :: t_log_g_lower_limit = t_log_L_upper_limit + 1 
      integer, parameter :: t_log_g_upper_limit = t_log_g_lower_limit + 1 
      integer, parameter :: t_power_nuc_burn_upper_limit = t_log_g_upper_limit + 1 
      integer, parameter :: t_power_h_burn_upper_limit = t_power_nuc_burn_upper_limit + 1 
      integer, parameter :: t_power_he_burn_upper_limit = t_power_h_burn_upper_limit + 1 
      integer, parameter :: t_power_z_burn_upper_limit = t_power_he_burn_upper_limit + 1 
      integer, parameter :: t_power_nuc_burn_lower_limit = t_power_z_burn_upper_limit + 1 
      integer, parameter :: t_power_h_burn_lower_limit = t_power_nuc_burn_lower_limit + 1 
      integer, parameter :: t_power_he_burn_lower_limit = t_power_h_burn_lower_limit + 1 
      integer, parameter :: t_power_z_burn_lower_limit = t_power_he_burn_lower_limit + 1 
      integer, parameter :: t_center_R_lower_limit = t_power_z_burn_lower_limit + 1 
      integer, parameter :: t_center_Ye_lower_limit = t_center_R_lower_limit + 1 
      integer, parameter :: t_fe_core_infall_limit = t_center_Ye_lower_limit + 1 
      integer, parameter :: t_non_fe_core_infall_limit = t_fe_core_infall_limit + 1 
      integer, parameter :: t_non_fe_core_rebound_limit = t_non_fe_core_infall_limit + 1 
      integer, parameter :: t_v_div_csound_max_limit = t_non_fe_core_rebound_limit + 1 
      integer, parameter :: t_v_div_csound_surf_limit = t_v_div_csound_max_limit + 1 
      integer, parameter :: t_gamma1_limit = t_v_div_csound_surf_limit + 1 
      integer, parameter :: t_Pgas_div_P_limit = t_gamma1_limit + 1 
      integer, parameter :: t_Lnuc_div_L_lower_limit = t_Pgas_div_P_limit + 1 
      integer, parameter :: t_Lnuc_div_L_upper_limit = t_Lnuc_div_L_lower_limit + 1 
      integer, parameter :: t_v_surf_div_v_kh_lower_limit = t_Lnuc_div_L_upper_limit + 1 
      integer, parameter :: t_v_surf_div_v_kh_upper_limit = t_v_surf_div_v_kh_lower_limit + 1 
      integer, parameter :: t_v_surf_div_v_esc_limit = t_v_surf_div_v_kh_upper_limit + 1 
      integer, parameter :: t_v_surf_kms_limit = t_v_surf_div_v_esc_limit + 1 
      integer, parameter :: t_Lnuc_div_L_zams_limit = t_v_surf_kms_limit + 1 
      integer, parameter :: t_phase_PreMS = t_Lnuc_div_L_zams_limit + 1 
      integer, parameter :: t_phase_ZAMS = t_phase_PreMS + 1 
      integer, parameter :: t_phase_IAMS = t_phase_ZAMS + 1 
      integer, parameter :: t_phase_TAMS = t_phase_IAMS + 1 
      integer, parameter :: t_phase_He_Burn = t_phase_TAMS + 1 
      integer, parameter :: t_phase_ZACHeB = t_phase_He_Burn + 1 
      integer, parameter :: t_phase_TACHeB = t_phase_ZACHeB + 1 
      integer, parameter :: t_phase_TP_AGB = t_phase_TACHeB + 1 
      integer, parameter :: t_phase_C_Burn = t_phase_TP_AGB + 1 
      integer, parameter :: t_phase_Ne_Burn = t_phase_C_Burn + 1 
      integer, parameter :: t_phase_O_Burn = t_phase_Ne_Burn + 1 
      integer, parameter :: t_phase_Si_Burn = t_phase_O_Burn + 1 
      integer, parameter :: t_phase_WDCS = t_phase_Si_Burn + 1       
      integer, parameter :: t_xa_central_lower_limit = t_phase_WDCS + 1 
      integer, parameter :: t_xa_central_upper_limit = t_xa_central_lower_limit + 1 
      integer, parameter :: t_xa_surface_lower_limit = t_xa_central_upper_limit + 1 
      integer, parameter :: t_xa_surface_upper_limit = t_xa_surface_lower_limit + 1 
      integer, parameter :: t_xa_average_lower_limit = t_xa_surface_upper_limit + 1 
      integer, parameter :: t_xa_average_upper_limit = t_xa_average_lower_limit + 1 
      integer, parameter :: t_surface_accel_div_grav_limit = t_xa_average_upper_limit + 1 
      integer, parameter :: t_adjust_mesh_failed = t_surface_accel_div_grav_limit + 1
      integer, parameter :: t_dt_is_zero = t_adjust_mesh_failed + 1 
      integer, parameter :: t_min_timestep_limit = t_dt_is_zero + 1 
      integer, parameter :: t_failed_prepare_for_new_try = t_min_timestep_limit + 1 
      integer, parameter :: t_negative_total_angular_momentum = t_failed_prepare_for_new_try + 1 
      integer, parameter :: t_max_number_retries = t_negative_total_angular_momentum + 1 
      integer, parameter :: t_redo_limit = t_max_number_retries + 1 
      integer, parameter :: t_solve_burn = t_redo_limit + 1 
      integer, parameter :: t_solve_hydro = t_solve_burn + 1 
      integer, parameter :: t_solve_mix = t_solve_hydro + 1 
      integer, parameter :: t_solve_omega_mix = t_solve_mix + 1 
      integer, parameter :: t_timestep_controller = t_solve_omega_mix + 1 
      integer, parameter :: t_relax_finished_okay = t_timestep_controller + 1
      integer, parameter :: t_cumulative_extra_heating_limit = t_relax_finished_okay + 1
      integer, parameter :: t_delta_total_energy = t_cumulative_extra_heating_limit + 1
      integer, parameter :: t_max_explicit_hydro_nsteps = t_delta_total_energy + 1
      integer, parameter :: t_center_entropy_limit = t_max_explicit_hydro_nsteps + 1
      integer, parameter :: t_center_entropy_lower_limit = t_center_entropy_limit + 1
      integer, parameter :: t_max_entropy_limit = t_center_entropy_lower_limit + 1
      integer, parameter :: t_max_entropy_lower_limit = t_max_entropy_limit + 1
      integer, parameter :: t_max_period_number = t_max_entropy_lower_limit + 1
      integer, parameter :: t_max_abs_rel_run_E_err = t_max_period_number + 1

      integer, parameter :: t_extras_check_model = t_max_abs_rel_run_E_err + 1
      integer, parameter :: t_extras_finish_step = t_extras_check_model + 1

      integer, parameter :: t_xtra1 = t_extras_finish_step + 1
      integer, parameter :: t_xtra2 = t_xtra1 + 1
      integer, parameter :: t_xtra3 = t_xtra2 + 1
      integer, parameter :: t_xtra4 = t_xtra3 + 1
      integer, parameter :: t_xtra5 = t_xtra4 + 1
      integer, parameter :: t_xtra6 = t_xtra5 + 1
      integer, parameter :: t_xtra7 = t_xtra6 + 1
      integer, parameter :: t_xtra8 = t_xtra7 + 1
      integer, parameter :: t_xtra9 = t_xtra8 + 1

      integer, parameter :: num_termination_codes = t_xtra9
      
      character (len=128) :: termination_code_str(num_termination_codes)
      
      
      ! auto_diff constants for solver variables
      ! used to access auto_diff_real_star_order1 d1Array
      integer, parameter :: i_lnd_m1 = 1
      integer, parameter :: i_lnd_00 = 2
      integer, parameter :: i_lnd_p1 = 3
      integer, parameter :: i_lnT_m1 = 4
      integer, parameter :: i_lnT_00 = 5
      integer, parameter :: i_lnT_p1 = 6
      integer, parameter :: i_w_m1 = 7
      integer, parameter :: i_w_00 = 8
      integer, parameter :: i_w_p1 = 9
      integer, parameter :: i_lnR_m1 = 10
      integer, parameter :: i_lnR_00 = 11
      integer, parameter :: i_lnR_p1 = 12
      integer, parameter :: i_v_m1 = 13
      integer, parameter :: i_v_00 = 14
      integer, parameter :: i_v_p1 = 15
      integer, parameter :: i_L_m1 = 16
      integer, parameter :: i_L_00 = 17
      integer, parameter :: i_L_p1 = 18
      integer, parameter :: i_hp_m1 = 19
      integer, parameter :: i_hp_00 = 20
      integer, parameter :: i_hp_p1 = 21
      integer, parameter :: i_w_div_wc_m1 = 22
      integer, parameter :: i_w_div_wc_00 = 23
      integer, parameter :: i_w_div_wc_p1 = 24
      integer, parameter :: i_jrot_m1 = 25
      integer, parameter :: i_jrot_00 = 26
      integer, parameter :: i_jrot_p1 = 27
      integer, parameter :: i_xtra1_m1 = 28
      integer, parameter :: i_xtra1_00 = 29
      integer, parameter :: i_xtra1_p1 = 30
      integer, parameter :: i_xtra2_m1 = 31
      integer, parameter :: i_xtra2_00 = 32
      integer, parameter :: i_xtra2_p1 = 33
      integer, parameter :: auto_diff_star_num_vars = i_xtra2_p1
      character (len=128) :: auto_diff_star_d1_names(auto_diff_star_num_vars)
      

      ! mesh cell types
      integer, parameter :: split_type = -1
      integer, parameter :: unchanged_type = 0
      integer, parameter :: merged_type = 1
      integer, parameter :: revised_type = 2
      
      
      ! solver operations
      integer, parameter :: struct_hydro_bit_pos = 0
      integer, parameter :: struct_thermo_bit_pos = 1
      integer, parameter :: burn_bit_pos = 2
      integer, parameter :: mix_bit_pos = 3
      
      ! relax rotation options
      integer, parameter :: relax_to_new_omega = 0
      integer, parameter :: relax_to_new_omega_div_omega_crit = 1
      integer, parameter :: relax_to_new_surface_rotation_v = 2

      
      ! values for result_reason
      integer, parameter :: result_reason_normal = 1
      integer, parameter :: dt_is_zero = 2 
         ! indicates that t+dt == t, probably because of round-off with tiny dt << t*1d-16.
      integer, parameter :: nonzero_ierr = 3 ! some routine returned with ierr /= 0
      integer, parameter :: hydro_failed_to_converge = 4
      integer, parameter :: do_burn_failed = 5
      integer, parameter :: diffusion_failed = 6
      integer, parameter :: too_many_steps_for_burn = 7
      integer, parameter :: too_many_steps_for_diffusion = 8
      integer, parameter :: too_many_steps_for_hydro = 9
      integer, parameter :: adjust_mesh_failed = 10
      integer, parameter :: adjust_mass_failed = 11
      integer, parameter :: core_dump_model_number = 12
      integer, parameter :: timestep_limits = 13
      integer, parameter :: variable_change_limits = 14
      integer, parameter :: explicit_hydro_failed = 15
      integer, parameter :: abs_rel_run_E_err = 16
      integer, parameter :: forced_stop = 17
      
      integer, parameter :: num_reasons = 17
      
      character (len=strlen) :: result_reason_str(num_reasons)
         
         
      character (len=strlen) :: rate_tables_dir_for_star, rates_cache_suffix_for_star

      character (len=strlen) :: compiler_name, compiler_version_name
      character (len=strlen) :: mesasdk_version_name, date


      ! debugging storage
      integer, parameter :: max_ndbg = 9
      
      ! bcyclic data
      type ulstore
         integer :: ul_size    ! size of umat1 & lmat1 (0 if not allocated)
         real(dp), pointer :: umat1(:), lmat1(:)
      end type ulstore
      
      type ulstore_qp
         integer :: ul_size    ! size of umat1 & lmat1 (0 if not allocated)
         real(qp), pointer :: umat1(:), lmat1(:)
      end type ulstore_qp

      integer, parameter :: max_generations = 3
      integer, parameter :: net_name_len = strlen
      integer, parameter :: name_len = 80           
      integer, parameter :: max_num_mixing_regions = 100
      integer, parameter :: maxlen_history_column_name = 80
      integer, parameter :: maxlen_profile_column_name = maxlen_history_column_name

      integer, parameter :: max_num_rsp_LINA_modes = 100

      integer, parameter :: max_num_accretion_species = 250

      integer, parameter :: max_num_diffusion_classes = 100
      
      integer, parameter :: num_termination_code_strings = 9
      
      integer, parameter :: num_xa_function = 9           
      integer, parameter :: num_mesh_logX = 9
      
      integer, parameter :: num_xa_central_limits = 9
      integer, parameter :: num_xa_surface_limits = 9
      integer, parameter :: num_xa_average_limits = 9

      integer, parameter :: NUM_OVERSHOOT_PARAM_SETS = 16
      integer, parameter :: NUM_PREDICTIVE_PARAM_SETS = 16
      integer, parameter :: NUM_SDOS_PARAM_SETS = 16
      
      integer, parameter :: terminal_iounit = 6

      integer, parameter :: max_num_profile_extras = 100
      integer, parameter :: max_num_history_extras = 100

      integer, parameter :: num_x_ctrls = 100
      
      integer, parameter :: GYRE_MODEL_VERSION = 101
      integer, parameter :: GSM_MODEL_VERSION = 110

      ! interfaces for procedure pointers
      abstract interface
      
      
      
      
         integer function id_only_function_interface(id)
            integer, intent(in) :: id
         end function id_only_function_interface
      
         subroutine id_only_subroutine_interface(id)
            integer, intent(in) :: id
         end subroutine id_only_subroutine_interface
      
      
         subroutine extras_startup_interface(id,restart,ierr)
            integer, intent(in) :: id
            integer, intent(out) :: ierr
            logical,intent(in) :: restart    
         end subroutine extras_startup_interface
         
         integer function extras_check_model_interface(id)
            integer, intent(in) :: id
         end function extras_check_model_interface
         
         integer function extras_start_step_interface(id)
            integer, intent(in) :: id
         end function extras_start_step_interface
         
         integer function extras_finish_step_interface(id)
            integer, intent(in) :: id
         end function extras_finish_step_interface

         subroutine extras_after_evolve_interface(id, ierr)
            integer, intent(in) :: id
            integer, intent(out) :: ierr
         end subroutine extras_after_evolve_interface
      

         subroutine other_wind_interface(id, Lsurf, Msurf, Rsurf, Tsurf, X, Y, Z, w, ierr)
            use const_def, only: dp
            integer, intent(in) :: id
            real(dp), intent(in) :: Lsurf, Msurf, Rsurf, Tsurf, X, Y, Z ! surface values (cgs)
            real(dp), intent(out) :: w ! wind in units of Msun/year (value is >= 0)
            integer, intent(out) :: ierr
         end subroutine other_wind_interface
      

         subroutine other_accreting_state_interface(id, total_specific_energy, accretion_pressure, accretion_density, ierr)
            use const_def, only: dp
            integer, intent(in) :: id
            real(dp), intent(out) :: total_specific_energy, accretion_pressure, accretion_density ! erg/g, erg/cm^3, g/cm^3 respectively
            integer, intent(out) :: ierr
         end subroutine other_accreting_state_interface
            
         subroutine other_before_struct_burn_mix_interface(id, dt, res)
            use const_def, only: dp
            integer, intent(in) :: id
            real(dp), intent(in) :: dt
            integer, intent(out) :: res ! keep_going, redo, retry, terminate
         end subroutine other_before_struct_burn_mix_interface
      
      
         subroutine other_after_struct_burn_mix_interface(id, dt, res)
            use const_def, only: dp
            integer, intent(in) :: id
            real(dp), intent(in) :: dt
            integer, intent(out) :: res ! keep_going, redo, retry, terminate
         end subroutine other_after_struct_burn_mix_interface
      
      
         subroutine other_adjust_mdot_interface(id, ierr)
            integer, intent(in) :: id
            integer, intent(out) :: ierr
         end subroutine other_adjust_mdot_interface
      
      
         subroutine other_after_solver_setmatrix_interface(id, ierr)
            integer, intent(in) :: id
            integer, intent(out) :: ierr
         end subroutine other_after_solver_setmatrix_interface
      
      
         subroutine other_j_for_adjust_J_lost_interface(id, starting_j_rot_surf, j_for_mass_loss, ierr)
            use const_def, only: dp
            integer, intent(in) :: id
            real(dp), intent(in) :: starting_j_rot_surf
            real(dp), intent(out) :: j_for_mass_loss
            integer, intent(out) :: ierr
         end subroutine other_j_for_adjust_J_lost_interface
      
      
         subroutine other_D_mix_interface(id, ierr)
            integer, intent(in) :: id
            integer, intent(out) :: ierr
         end subroutine other_D_mix_interface
      
      
         subroutine other_am_mixing_interface(id, ierr)
            integer, intent(in) :: id
            integer, intent(out) :: ierr
         end subroutine other_am_mixing_interface


         subroutine other_diffusion_coefficients_interface( &
               id, k, nc, m, rho, T, A, X, Z, C, charge, na, &
               Ddiff, Kdiff, Zdiff, Zdiff1, Zdiff2, Ath)
            use const_def, only: dp
            integer, intent(in) :: id, k, nc, m  
            real(dp), intent(in) :: rho, T, charge(m), na(m)
            real(dp), intent(in), dimension(:) :: A, X, Z, C ! (m)
            real(dp), intent(inout), dimension(m,m) :: &
               Ddiff, Kdiff, Zdiff, Zdiff1, Zdiff2, Ath
         end subroutine other_diffusion_coefficients_interface
      
         subroutine other_remove_surface_interface(id, ierr, k)
            integer, intent(in) :: id
            integer, intent(out) :: ierr
            integer, intent(out) :: k
         end subroutine other_remove_surface_interface
      
         subroutine other_momentum_interface(id, ierr)
            integer, intent(in) :: id
            integer, intent(out) :: ierr
         end subroutine other_momentum_interface      
      
         subroutine other_momentum_implicit_interface(id, ierr)
            integer, intent(in) :: id
            integer, intent(out) :: ierr
         end subroutine other_momentum_implicit_interface
      
      
         subroutine other_pressure_interface(id, ierr)
            integer, intent(in) :: id
            integer, intent(out) :: ierr
         end subroutine other_pressure_interface
      
      
         subroutine other_energy_interface(id, ierr)
            integer, intent(in) :: id
            integer, intent(out) :: ierr
         end subroutine other_energy_interface
      
      
         subroutine other_energy_implicit_interface(id, ierr)
            integer, intent(in) :: id
            integer, intent(out) :: ierr
         end subroutine other_energy_implicit_interface
      
      
         subroutine other_brunt_interface(id, ierr)
            integer, intent(in) :: id
            integer, intent(out) :: ierr
         end subroutine other_brunt_interface
      
      
         subroutine other_brunt_smoothing_interface(id, ierr)
            integer, intent(in) :: id
            integer, intent(out) :: ierr
         end subroutine other_brunt_smoothing_interface
      
      
         subroutine other_build_initial_model_interface(id, ierr)
            integer, intent(in) :: id
            integer, intent(out) :: ierr
         end subroutine other_build_initial_model_interface
      
         subroutine other_eval_fp_ft_interface( &
               id, nz, xm, r, rho, aw, ft, fp, r_polar, r_equatorial, report_ierr, ierr)
            use const_def, only: dp
            use auto_diff
            integer, intent(in) :: id
            integer, intent(in) :: nz
            real(dp), intent(in) :: aw(:), r(:), rho(:), xm(:)
            type(auto_diff_real_star_order1), intent(out) :: ft(:), fp(:)
            real(dp), intent(inout) :: r_polar(:), r_equatorial(:)
            logical, intent(in) :: report_ierr
            integer, intent(out) :: ierr
         end subroutine other_eval_fp_ft_interface
      
         subroutine other_eval_i_rot_interface(id,k,r00,w_div_w_crit_roche, i_rot)
            use const_def, only: dp
            use auto_diff
            integer, intent(in) :: id,k
            real(dp), intent(in) :: r00,w_div_w_crit_roche
            type(auto_diff_real_star_order1), intent(out) :: i_rot
         end subroutine other_eval_i_rot_interface
      
         subroutine other_torque_interface(id, ierr)
            integer, intent(in) :: id
            integer, intent(out) :: ierr
         end subroutine other_torque_interface
      
      
         subroutine other_torque_implicit_interface(id, ierr)
            integer, intent(in) :: id
            integer, intent(out) :: ierr
         end subroutine other_torque_implicit_interface
      
      
         subroutine other_cgrav_interface(id, ierr)
            integer, intent(in) :: id
            integer, intent(out) :: ierr
         end subroutine other_cgrav_interface
      
      
         subroutine other_overshooting_scheme_interface(id, i, j, k_a, k_b, D, vc, ierr)
            use const_def, only: dp
            integer, intent(in) :: id, i, j
            integer, intent(out) :: k_a, k_b
            real(dp), intent(out), dimension(:) :: D, vc
            integer, intent(out) :: ierr
         end subroutine other_overshooting_scheme_interface
      
      
         subroutine other_mesh_delta_coeff_factor_interface(id, ierr)
            use const_def, only: dp
            integer, intent(in) :: id
            integer, intent(out) :: ierr
         end subroutine other_mesh_delta_coeff_factor_interface
      
      
         subroutine other_alpha_mlt_interface(id, ierr)
            integer, intent(in) :: id
            integer, intent(out) :: ierr
         end subroutine other_alpha_mlt_interface
      
      
         subroutine other_eps_grav_interface(id, k, dt, ierr)
            use const_def, only: dp
            integer, intent(in) :: id, k
            real(dp), intent(in) :: dt
            integer, intent(out) :: ierr
         end subroutine other_eps_grav_interface
      
      
         subroutine other_rsp_build_model_interface(id, ierr)
            integer, intent(in) :: id
            integer, intent(out) :: ierr
         end subroutine other_rsp_build_model_interface
      
      
         subroutine other_rsp_linear_analysis_interface(id, restart, ierr)
            integer, intent(in) :: id
            logical, intent(in) :: restart
            integer, intent(out) :: ierr
         end subroutine other_rsp_linear_analysis_interface
      
      
         subroutine other_gradr_factor_interface(id, ierr)
            integer, intent(in) :: id
            integer, intent(out) :: ierr
         end subroutine other_gradr_factor_interface
      
      
         subroutine other_opacity_factor_interface(id, ierr)
            integer, intent(in) :: id
            integer, intent(out) :: ierr
         end subroutine other_opacity_factor_interface


         subroutine other_diffusion_factor_interface(id, ierr)
            integer, intent(in) :: id
            integer, intent(out) :: ierr
         end subroutine other_diffusion_factor_interface

      
         subroutine other_solver_monitor_interface( &
               id, iter, passed_tol_tests, &
               correction_norm, max_correction, &
               residual_norm, max_residual, ierr)
            use const_def, only: dp
            integer, intent(in) :: id, iter 
            logical, intent(in) :: passed_tol_tests
            real(dp), intent(in) :: correction_norm, max_correction, &
               residual_norm, max_residual
            integer, intent(out) :: ierr
         end subroutine other_solver_monitor_interface


         subroutine other_neu_interface(  &
               id, k, T, log10_T, Rho, log10_Rho, abar, zbar, log10_Tlim, flags, &
               loss, sources, ierr)
            use const_def, only: dp
            use neu_def, only: num_neu_types, num_neu_rvs
            integer, intent(in) :: id ! id for star         
            integer, intent(in) :: k ! cell number or 0 if not for a particular cell         
            real(dp), intent(in) :: T ! temperature
            real(dp), intent(in) :: log10_T ! log10 of temperature
            real(dp), intent(in) :: Rho ! density
            real(dp), intent(in) :: log10_Rho ! log10 of density
            real(dp), intent(in) :: abar ! mean atomic weight
            real(dp), intent(in) :: zbar ! mean charge
            real(dp), intent(in) :: log10_Tlim 
            logical, intent(inout) :: flags(num_neu_types) ! true if should include the type of loss
            real(dp), intent(inout) :: loss(num_neu_rvs) ! total from all sources
            real(dp), intent(inout) :: sources(num_neu_types, num_neu_rvs)
            integer, intent(out) :: ierr
         end subroutine other_neu_interface
         
         
         subroutine other_net_get_interface(  &
               id, k, net_handle, just_dxdt, n, num_isos, num_reactions,  &
               x, temp, log10temp, rho, log10rho,  &
               abar, zbar, z2bar, ye, eta, d_eta_dlnT, d_eta_dlnRho, &
               rate_factors, weak_rate_factor, &
               reaction_Qs, reaction_neuQs, &
               eps_nuc, d_eps_nuc_dRho, d_eps_nuc_dT, d_eps_nuc_dx,  &
               dxdt, d_dxdt_dRho, d_dxdt_dT, d_dxdt_dx,  &
               screening_mode,  &
               eps_nuc_categories, eps_neu_total, &
               lwork, work, ierr)      
            use const_def, only: dp
            use net_def, only: Net_Info
            integer, intent(in) :: id ! id for star         
            integer, intent(in) :: k ! cell number or 0 if not for a particular cell         
            integer, intent(in) :: net_handle
            logical, intent(in) :: just_dxdt
            type (Net_Info) :: n
            integer, intent(in) :: num_isos
            integer, intent(in) :: num_reactions
            real(dp), intent(in)  :: x(:) ! (num_isos)
            real(dp), intent(in)  :: temp, log10temp ! log10 of temp
            real(dp), intent(in)  :: rho, log10rho ! log10 of rho
            real(dp), intent(in)  :: abar  ! mean number of nucleons per nucleus
            real(dp), intent(in)  :: zbar  ! mean charge per nucleus
            real(dp), intent(in)  :: z2bar ! mean charge squared per nucleus
            real(dp), intent(in)  :: ye    
            real(dp), intent(in)  :: eta, d_eta_dlnT, d_eta_dlnRho ! electron degeneracy from eos.
            real(dp), intent(in), pointer :: rate_factors(:) ! (num_reactions)
            real(dp), intent(in) :: weak_rate_factor
            real(dp), pointer, intent(in) :: reaction_Qs(:) ! (rates_reaction_id_max)
            real(dp), pointer, intent(in) :: reaction_neuQs(:) ! (rates_reaction_id_max)
            real(dp), intent(out) :: eps_nuc ! ergs/g/s from burning including losses by reaction neutrinos
            real(dp), intent(out) :: d_eps_nuc_dT
            real(dp), intent(out) :: d_eps_nuc_dRho
            real(dp), intent(inout) :: d_eps_nuc_dx(:) ! (num_isos)       
            real(dp), intent(inout) :: dxdt(:) ! (num_isos)
            real(dp), intent(inout) :: d_dxdt_dRho(:) ! (num_isos)
            real(dp), intent(inout) :: d_dxdt_dT(:) ! (num_isos)
            real(dp), intent(inout) :: d_dxdt_dx(:,:) ! (num_isos, num_isos)            
            real(dp), intent(inout) :: eps_nuc_categories(:) ! (num_categories)
            real(dp), intent(out) :: eps_neu_total ! ergs/g/s neutrinos from weak reactions
            integer, intent(in) :: screening_mode
            integer, intent(in) :: lwork ! size of work >= result from calling net_work_size
            real(dp), pointer :: work(:) ! (lwork)
            integer, intent(out) :: ierr ! 0 means okay
         end subroutine other_net_get_interface


         subroutine other_split_burn_interface(  &
            id, k, net_handle, eos_handle, num_isos, num_reactions, t_start, t_end, starting_x, &
            num_times_for_interpolation, times, log10Ts_f1, log10Rhos_f1, etas_f1, &
            dxdt_source_term, rate_factors, &
            weak_rate_factor, reaction_Qs, reaction_neuQs, &
            screening_mode, &
            stptry, max_steps, eps, odescal, &
            use_pivoting, trace, dbg, burner_finish_substep, &
            burn_lwork, burn_work_array, &
            net_lwork, net_work_array, &
            ! results
            ending_x, eps_nuc_categories, avg_eps_nuc, eps_neu_total, &
            nfcn, njac, nstep, naccpt, nrejct, ierr)
            use net_def
            
            integer, intent(in) :: net_handle, eos_handle, id
            integer, intent(in) :: num_isos
            integer, intent(in) :: num_reactions
            real(dp), intent(in) :: t_start, t_end, starting_x(:) ! (num_isos)
            
            integer, intent(in) :: num_times_for_interpolation 
               ! ending time is times(num_times); starting time is 0
            real(dp), pointer, intent(in) :: times(:) ! (num_times) 
            real(dp), pointer, intent(in) :: log10Ts_f1(:) ! =(4,numtimes) interpolant for log10T(time)
            real(dp), pointer, intent(in) :: log10Rhos_f1(:) ! =(4,numtimes) interpolant for log10Rho(time)
            real(dp), pointer, intent(in) :: etas_f1(:) ! =(4,numtimes) interpolant for eta(time)
            real(dp), pointer, intent(in) :: dxdt_source_term(:) ! (num_isos)  or null if no source term.
            real(dp), intent(in), pointer :: rate_factors(:) ! (num_reactions)
            real(dp), intent(in) :: weak_rate_factor
            real(dp), pointer, intent(in) :: reaction_Qs(:) ! (rates_reaction_id_max)
            real(dp), pointer, intent(in) :: reaction_neuQs(:) ! (rates_reaction_id_max)
            integer, intent(in) :: screening_mode ! see screen_def
            real(dp), intent(in) :: stptry ! try this for 1st step.  0 means try in 1 step.
            integer, intent(in) :: max_steps ! maximal number of allowed steps.
            real(dp), intent(in) :: eps, odescal ! tolerances.  e.g., set both to 1d-6
            logical, intent(in) :: use_pivoting ! for matrix solves
            logical, intent(in) :: trace, dbg
            interface
               include 'burner_finish_substep.inc'
            end interface
            integer, intent(in) :: net_lwork, burn_lwork
            real(dp), intent(inout), pointer :: burn_work_array(:) ! (burn_lwork)
            real(dp), intent(inout), pointer :: net_work_array(:) ! (net_lwork)
   
            ! These should be set for the output
            real(dp), intent(inout) :: ending_x(:) ! (num_isos)
            real(dp), intent(inout) :: eps_nuc_categories(:) ! (num_categories)
            real(dp), intent(out) :: avg_eps_nuc, eps_neu_total
            integer, intent(out) :: nfcn    ! number of function evaluations
            integer, intent(out) :: njac    ! number of jacobian evaluations
            integer, intent(out) :: nstep   ! number of computed steps
            integer, intent(out) :: naccpt  ! number of accepted steps
            integer, intent(out) :: nrejct  ! number of rejected steps
            integer, intent(out) :: ierr
            integer,intent(in) :: k
      
         end subroutine other_split_burn_interface

         
         subroutine other_adjust_mlt_gradT_fraction_interface(id, ierr)
            integer, intent(in) :: id
            integer, intent(out) :: ierr
         end subroutine other_adjust_mlt_gradT_fraction_interface


         subroutine other_mlt_results_interface(id, k, MLT_option, &  ! NOTE: k=0 is a valid arg
               r, L, T, P, opacity, rho, chiRho, chiT, Cp, gradr, grada, scale_height, &
               iso, XH1, cgrav, m, gradL_composition_term, mixing_length_alpha, &
               alpha_semiconvection, thermohaline_coeff, &
               mixing_type, gradT, Y_face, conv_vel, D, Gamma, ierr)
            use const_def, only: dp
            use auto_diff
            integer, intent(in) :: id
            integer, intent(in) :: k
            character (len=*), intent(in) :: MLT_option
            type(auto_diff_real_star_order1), intent(in) :: &
               r, L, T, P, opacity, rho, chiRho, chiT, Cp, gradr, grada, scale_height
            integer, intent(in) :: iso
            real(dp), intent(in) :: &
               XH1, cgrav, m, gradL_composition_term, &
               mixing_length_alpha, alpha_semiconvection, thermohaline_coeff
            integer, intent(out) :: mixing_type
            type(auto_diff_real_star_order1), intent(out) :: &
               gradT, Y_face, conv_vel, D, Gamma
            integer, intent(out) :: ierr
         end subroutine other_mlt_results_interface

      
         subroutine other_after_set_mixing_info_interface(id, ierr)
            integer, intent(in) :: id
            integer, intent(out) :: ierr
         end subroutine other_after_set_mixing_info_interface
      
      
         subroutine other_diffusion_interface(id, dt, ierr)
            use const_def, only: dp
            integer, intent(in) :: id
            real(dp), intent(in) :: dt 
            integer, intent(out) :: ierr
         end subroutine other_diffusion_interface


         subroutine other_kap_get_interface( &
               id, k, handle, species, chem_id, net_iso, xa, &
               log10_rho, log10_T, &
               lnfree_e, d_lnfree_e_dlnRho, d_lnfree_e_dlnT, &
               eta, d_eta_dlnRho, d_eta_dlnT, &
               kap_fracs, kap, dln_kap_dlnRho, dln_kap_dlnT, dln_kap_dxa, ierr)
      
            use const_def, only: dp
            use kap_def, only: num_kap_fracs
 
            ! INPUT
            integer, intent(in) :: id ! star id if available; 0 otherwise
            integer, intent(in) :: k ! cell number or 0 if not for a particular cell         
            integer, intent(in) :: handle ! from alloc_kap_handle
            integer, intent(in) :: species
            integer, pointer :: chem_id(:) ! maps species to chem id
               ! index from 1 to species
               ! value is between 1 and num_chem_isos         
            integer, pointer :: net_iso(:) ! maps chem id to species number
               ! index from 1 to num_chem_isos (defined in chem_def)
               ! value is 0 if the iso is not in the current net
               ! else is value between 1 and number of species in current net
            real(dp), intent(in) :: xa(:) ! mass fractions
            real(dp), intent(in) :: log10_rho ! density
            real(dp), intent(in) :: log10_T ! temperature
            real(dp), intent(in) :: lnfree_e, d_lnfree_e_dlnRho, d_lnfree_e_dlnT
               ! free_e := total combined number per nucleon of free electrons and positrons
            real(dp), intent(in) :: eta, d_eta_dlnRho, d_eta_dlnT
               ! eta := electron degeneracy parameter
         
            ! OUTPUT
            real(dp), intent(out) :: kap_fracs(num_kap_fracs)
            real(dp), intent(out) :: kap ! opacity
            real(dp), intent(out) :: dln_kap_dlnRho ! partial derivative at constant T
            real(dp), intent(out) :: dln_kap_dlnT   ! partial derivative at constant Rho
            real(dp), intent(out) :: dln_kap_dxa(:) ! partial derivative w.r.t. to species
            integer, intent(out) :: ierr ! 0 means AOK.

         end subroutine other_kap_get_interface
      
      
         subroutine other_kap_get_op_mono_interface( &
               handle, zbar, log10_rho, log10_T, &
               ! args for op_mono
               use_op_mono_alt_get_kap, &
               nel, izzp, fap, fac, screening, umesh, semesh, ff, rs, &
               ! output
               kap, dlnkap_dlnRho, dlnkap_dlnT, ierr)
            use const_def, only: dp
            integer, intent(in) :: handle ! from alloc_kap_handle
            real(dp), intent(in) :: zbar
            real(dp), intent(in) :: log10_rho ! the density
            real(dp), intent(in) :: log10_T ! the temperature
            ! args for op_mono_get_kap
            logical, intent(in) :: use_op_mono_alt_get_kap
            integer, intent(in) :: nel
            integer, intent(in) :: izzp(:) ! (nel)
            real(dp), intent(in) :: fap(:) ! (nel) number fractions of elements
            real(dp), intent(in) :: fac(:) ! (nel) scale factors for element opacity
            logical, intent(in) :: screening
            ! work arrays
            real, pointer :: umesh(:), semesh(:), ff(:,:,:,:), rs(:,:,:)
               ! umesh(nptot)
               ! semesh(nptot)
               ! ff(nptot, ipe, 4, 4)
               ! rs(nptot, 4, 4)
            ! output
            real(dp), intent(out) :: kap ! opacity
            real(dp), intent(out) :: dlnkap_dlnRho
            real(dp), intent(out) :: dlnkap_dlnT
            integer, intent(out) :: ierr ! 0 means AOK.
         end subroutine other_kap_get_op_mono_interface
      
      
         subroutine other_pgstar_plots_info_interface(id, ierr)
            integer, intent(in) :: id
            integer, intent(out) :: ierr
         end subroutine other_pgstar_plots_info_interface

      
         subroutine how_many_other_mesh_fcns_interface(id, n)
            integer, intent(in) :: id
            integer, intent(out) :: n
         end subroutine how_many_other_mesh_fcns_interface
      
      
         subroutine other_mesh_fcn_data_interface( &
               id, nfcns, names, gval_is_xa_function, vals1, ierr)
            use const_def, only: dp
            integer, intent(in) :: id
            integer, intent(in) :: nfcns
            character (len=*) :: names(:)
            logical, intent(out) :: gval_is_xa_function(:) ! (nfcns)
            real(dp), pointer :: vals1(:) ! =(nz, nfcns)
            integer, intent(out) :: ierr
         end subroutine other_mesh_fcn_data_interface
      
      
         subroutine other_write_pulsation_info_interface( &
               id, add_center_point, keep_surface_point, &
               add_atmosphere, filename, ierr)
            integer, intent(in) :: id
            logical, intent(in) :: &
               add_center_point, keep_surface_point, &
               add_atmosphere
            character (len=*), intent(in) :: filename
            integer, intent(out) :: ierr 
         end subroutine other_write_pulsation_info_interface
      
      
         subroutine other_get_pulsation_info_interface(id, &
               add_center_point, keep_surface_point, &
               add_atmosphere, nn, iconst, ivar, glob, var, ierr)
            use const_def, only: dp
            integer, intent(in) :: id
            logical, intent(in) :: &
               add_center_point, keep_surface_point, &
               add_atmosphere
            integer, intent(out) :: nn, iconst, ivar
            real(dp), pointer :: glob(:) ! (iconst) -- will be allocated
            real(dp), pointer :: var(:,:) ! (ivar,nn) -- will be allocated
            integer, intent(out) :: ierr
         end subroutine other_get_pulsation_info_interface
      
      
         subroutine other_surface_PT_interface(id, &
               skip_partials, &
               lnT_surf, dlnT_dL, dlnT_dlnR, dlnT_dlnM, dlnT_dlnkap, &
               lnP_surf, dlnP_dL, dlnP_dlnR, dlnP_dlnM, dlnP_dlnkap, ierr)
            use const_def, only: dp
            integer, intent(in) :: id
            logical, intent(in) :: skip_partials
            real(dp), intent(out) :: &
               lnT_surf, dlnT_dL, dlnT_dlnR, dlnT_dlnM, dlnT_dlnkap, &
               lnP_surf, dlnP_dL, dlnP_dlnR, dlnP_dlnM, dlnP_dlnkap
            integer, intent(out) :: ierr
         end subroutine other_surface_PT_interface
      
      
         subroutine other_photo_write_interface(id, iounit)
            integer, intent(in) :: id, iounit
         end subroutine other_photo_write_interface
      
      
         subroutine other_photo_read_interface(id, iounit, ierr)
            integer, intent(in) :: id, iounit
            integer, intent(out) :: ierr
         end subroutine other_photo_read_interface
      
      
         subroutine other_new_generation_interface(id, ierr)
            integer, intent(in) :: id
            integer, intent(out) :: ierr
         end subroutine other_new_generation_interface
      
      
         subroutine other_set_current_to_old_interface(id)
            integer, intent(in) :: id
         end subroutine other_set_current_to_old_interface
      
      
         subroutine other_set_pgstar_controls_interface(id)
            integer, intent(in) :: id
         end subroutine other_set_pgstar_controls_interface
      
      
         subroutine other_astero_freq_corr_interface(id, ierr)
            integer, intent(in) :: id
            integer, intent(out) :: ierr
         end subroutine other_astero_freq_corr_interface
      
      
         integer function other_timestep_limit_interface( &
            id, skip_hard_limit, dt, dt_limit_ratio)
            ! return keep_going, retry, or terminate
            use const_def, only: dp
            integer, intent(in) :: id
            logical, intent(in) :: skip_hard_limit
            real(dp), intent(in) :: dt
            real(dp), intent(inout) :: dt_limit_ratio
         end function other_timestep_limit_interface

      
         integer function how_many_extra_header_items_interface(id)
            integer, intent(in) :: id
         end function how_many_extra_header_items_interface


         subroutine data_for_extra_header_items_interface( &
               id, n, names, vals, ierr)
            use const_def, only: dp
            integer, intent(in) :: id, n
            character (len=80) :: names(n)
            real(dp) :: vals(n)
            integer, intent(out) :: ierr
         end subroutine data_for_extra_header_items_interface


         integer function how_many_extra_history_columns_interface(id)
            integer, intent(in) :: id
         end function how_many_extra_history_columns_interface
      
      
         subroutine data_for_extra_history_columns_interface( &
               id, n, names, vals, ierr)
            use const_def, only: dp
            integer, intent(in) :: id, n
            character (len=80) :: names(n)
            real(dp) :: vals(n)
            integer, intent(out) :: ierr
         end subroutine data_for_extra_history_columns_interface

      
         integer function how_many_extra_profile_columns_interface(id)
            integer, intent(in) :: id
         end function how_many_extra_profile_columns_interface
      
      
         subroutine data_for_extra_profile_columns_interface( &
               id, n, nz, names, vals, ierr)
            use const_def, only: dp
            integer, intent(in) :: id, n, nz
            character (len=80) :: names(n)
            real(dp) :: vals(nz,n)
            integer, intent(out) :: ierr
         end subroutine data_for_extra_profile_columns_interface

      
         integer function how_many_binary_history_columns_interface(binary_id)
            integer, intent(in) :: binary_id
         end function how_many_binary_history_columns_interface
      
      
         subroutine data_for_binary_history_columns_interface( &
               binary_id, n, names, vals, ierr)
            use const_def, only: dp
            integer, intent(in) :: binary_id, n
            character (len=80) :: names(n)
            real(dp) :: vals(n)
            integer, intent(out) :: ierr
         end subroutine data_for_binary_history_columns_interface
           

         integer function how_many_extra_binary_history_columns_interface(binary_id)
            integer, intent(in) :: binary_id
         end function how_many_extra_binary_history_columns_interface


         subroutine data_for_extra_binary_history_columns_interface( &
               binary_id, n, names, vals, ierr)
            use const_def, only: dp
            integer, intent(in) :: binary_id, n
            character (len=80) :: names(n)
            real(dp) :: vals(n)
            integer, intent(out) :: ierr
         end subroutine data_for_extra_binary_history_columns_interface


         integer function how_many_extra_binary_history_header_items_interface(binary_id)
            integer, intent(in) :: binary_id
         end function how_many_extra_binary_history_header_items_interface


         subroutine data_for_extra_binary_history_header_items_interface( &
               binary_id, n, names, vals, ierr)
            use const_def, only: dp
            integer, intent(in) :: binary_id, n
            character (len=80) :: names(n)
            real(dp) :: vals(n)
            integer, intent(out) :: ierr
         end subroutine data_for_extra_binary_history_header_items_interface


         subroutine other_export_pulse_data_interface (id, data_format, filename, &
              add_atmosphere, keep_surface_point, add_center_point, ierr)
           integer, intent(in)      :: id
           character(*), intent(in) :: data_format
           character(*), intent(in) :: filename
           logical, intent(in)      :: add_atmosphere
           logical, intent(in)      :: keep_surface_point
           logical, intent(in)      :: add_center_point
           integer, intent(out)     :: ierr
         end subroutine other_export_pulse_data_interface

         
         subroutine other_get_pulse_data_interface (id, data_format, &
              add_atmosphere, keep_surface_point, add_center_point, global_data, point_data, ierr)
           use const_def, only: dp
           integer, intent(in)                :: id
           character(*), intent(in)           :: data_format
           logical, intent(in)                :: add_atmosphere
           logical, intent(in)                :: keep_surface_point
           logical, intent(in)                :: add_center_point
           real(dp), allocatable, intent(out) :: global_data(:)
           real(dp), allocatable, intent(out) :: point_data(:,:)
           integer, intent(out)               :: ierr
         end subroutine other_get_pulse_data_interface

         
         subroutine other_edit_pulse_data_interface (id, data_format, global_data, point_data, ierr)
           use const_def, only: dp
           integer, intent(in)      :: id
           character(*), intent(in) :: data_format
           real(dp), intent(inout)  :: global_data(:)
           real(dp), intent(inout)  :: point_data(:,:)
           integer, intent(out)     :: ierr
         end subroutine other_edit_pulse_data_interface
      
      
         subroutine binary_other_torque_implicit_interface(id, ierr)
            integer, intent(in) :: id
            integer, intent(out) :: ierr
         end subroutine binary_other_torque_implicit_interface
      
      
         subroutine binary_other_torque_interface(id, ierr)
            integer, intent(in) :: id
            integer, intent(out) :: ierr
         end subroutine binary_other_torque_interface
                  
      end interface

