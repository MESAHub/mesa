! specifications for starting model
real(dp) :: m1
real(dp) :: m2
real(dp) :: initial_period_in_days
real(dp) :: initial_separation_in_Rsuns
real(dp) :: initial_eccentricity

! controls for output
character (len = strlen) :: history_name
integer :: history_interval
logical :: append_to_star_history
character (len = strlen) :: log_directory
character (len = strlen) :: history_dbl_format
character (len = strlen) :: history_int_format
character (len = strlen) :: history_txt_format
integer :: photo_interval
integer :: photo_digits
character (len = strlen) :: photo_directory
integer :: terminal_interval
integer :: write_header_frequency
character (len = strlen) :: extra_binary_terminal_output_file

! timestep controls
real(dp) :: fm
real(dp) :: fm_hard
real(dp) :: fa
real(dp) :: fa_hard
real(dp) :: fr
real(dp) :: fr_hard
real(dp) :: fj
real(dp) :: fj_hard
real(dp) :: fe
real(dp) :: fe_hard
real(dp) :: fm_limit
real(dp) :: fr_limit
real(dp) :: fe_limit
real(dp) :: fr_dt_limit
real(dp) :: fdm
real(dp) :: fdm_hard
real(dp) :: dt_softening_factor
real(dp) :: varcontrol_case_a
real(dp) :: varcontrol_case_b
real(dp) :: varcontrol_ms
real(dp) :: varcontrol_post_ms
real(dp) :: dt_reduction_factor_for_j

! when to stop
real(dp) :: accretor_overflow_terminate
logical :: terminate_if_initial_overflow
logical :: terminate_if_L2_overflow


! mass transfer controls
real(dp) :: mass_transfer_alpha
real(dp) :: mass_transfer_beta
real(dp) :: mass_transfer_delta
real(dp) :: mass_transfer_gamma
logical :: limit_retention_by_mdot_edd
logical :: use_es_opacity_for_mdot_edd
real(dp) :: use_this_for_mdot_edd_eta
logical :: use_radiation_corrected_transfer_rate
real(dp) :: initial_bh_spin
real(dp) :: use_this_for_mdot_edd
character (len = strlen) :: mdot_scheme
real(dp) :: cur_mdot_frac
real(dp) :: max_explicit_abs_mdot
real(dp) :: max_implicit_abs_mdot
logical :: report_rlo_solver_progress
integer :: max_tries_to_achieve
character (len = strlen) :: solver_type
real(dp) :: implicit_scheme_tolerance
real(dp) :: implicit_scheme_tiny_factor
real(dp) :: initial_change_factor
real(dp) :: change_factor_fraction
real(dp) :: implicit_lambda
real(dp) :: max_change_factor
real(dp) :: min_change_factor
integer :: num_tries_for_increase_change_factor
real(dp) :: change_factor_increase
real(dp) :: starting_mdot
real(dp) :: roche_min_mdot
real(dp) :: min_mdot_for_implicit
logical :: do_enhance_wind_1, do_enhance_wind_2
real(dp) :: tout_B_wind_1, tout_B_wind_2
logical :: do_wind_mass_transfer_1, do_wind_mass_transfer_2
real(dp) :: wind_BH_alpha_1, wind_BH_alpha_2
real(dp) :: max_wind_transfer_fraction_1, max_wind_transfer_fraction_2

! orbital jdot controls
logical :: do_jdot_gr
logical :: do_jdot_ml
logical :: do_jdot_ls
logical :: do_jdot_missing_wind
logical :: do_jdot_mb
logical :: include_accretor_mb
real(dp) :: magnetic_braking_gamma
logical :: keep_mb_on
real(dp) :: jdot_mb_min_qconv_env, jdot_mb_max_qconv_env, jdot_mb_max_qconv_core, &
      jdot_mb_qlim_for_check_rad_core, jdot_mb_qlim_for_check_conv_env
logical :: jdot_mb_scale_for_low_qconv_env
real(dp) :: jdot_mb_mass_frac_for_scale
real(dp) :: jdot_multiplier

! rotation and sync controls
logical :: do_j_accretion
logical :: do_tidal_sync
character (len = strlen) :: sync_type_1, sync_type_2
character (len = strlen) :: sync_mode_1, sync_mode_2
real(dp) :: Ftid_1, Ftid_2
logical :: do_initial_orbit_sync_1, do_initial_orbit_sync_2
real(dp) :: tidal_reduction

! eccentricity controls
logical :: do_tidal_circ
character (len = strlen) :: circ_type_1, circ_type_2
logical :: use_eccentricity_enhancement
real(dp) :: max_abs_edot_tidal
real(dp) :: max_abs_edot_enhance
real(dp) :: min_eccentricity
real(dp) :: max_eccentricity
integer :: anomaly_steps

! irradiation controls
logical :: accretion_powered_irradiation
real(dp) :: col_depth_for_eps_extra
logical :: use_accretor_luminosity_for_irrad
real(dp) :: irrad_flux_at_std_distance
real(dp) :: std_distance_for_irradiation
real(dp) :: max_F_irr

! common envelope controls
real(dp) :: CE_alpha
real(dp) :: CE_alpha_th
real(dp) :: CE_alpha_core
real(dp) :: CE_mass_loss_rate_high
real(dp) :: CE_mass_loss_rate_low
real(dp) :: CE_rel_rlo_for_detachment
real(dp) :: CE_years_detached_to_terminate
logical :: CE_begin_at_max_implicit_abs_mdot
real(dp) :: CE_xa_diff_to_terminate
logical :: CE_terminate_when_core_overflows
real(dp) :: CE_min_period_in_minutes
real(dp) :: CE_energy_factor_HII_toHI
real(dp) :: CE_energy_factor_HeII_toHeI
real(dp) :: CE_energy_factor_HeIII_toHeII
real(dp) :: CE_energy_factor_H2
real(dp) :: CE_fixed_lambda


! miscellaneous controls
logical :: keep_donor_fixed
real(dp) :: mdot_limit_donor_switch
logical :: use_other_rlo_mdot
logical :: use_other_check_implicit_rlo
logical :: use_other_implicit_function_to_solve
logical :: use_other_tsync
logical :: use_other_sync_spin_to_orbit
logical :: use_other_mdot_edd
logical :: use_other_adjust_mdots
logical :: use_other_accreted_material_j
logical :: use_other_jdot_gr
logical :: use_other_jdot_ml
logical :: use_other_jdot_ls
logical :: use_other_jdot_missing_wind
logical :: use_other_jdot_mb
logical :: use_other_extra_jdot
logical :: use_other_binary_wind_transfer
logical :: use_other_edot_tidal
logical :: use_other_edot_enhance
logical :: use_other_extra_edot
logical :: use_other_CE_init
logical :: use_other_CE_rlo_mdot
logical :: use_other_CE_binary_evolve_step
logical :: use_other_CE_binary_finish_step
logical :: use_other_e2
logical :: use_other_pgbinary_plots

! end of controls
