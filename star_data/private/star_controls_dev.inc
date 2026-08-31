         
         logical :: compare_TDC_to_MLT
         logical :: use_TDC_Y_face_seeded_newton
         character(32) :: hydro_matrix_solver

         logical :: TDC_adjust_mass_fallback_to_mlt
         logical :: TDC_use_dynamical_gradL
         logical :: include_mlt_Pturb_in_thermodynamic_gradients
         logical :: use_TDC_enthalpy_flux_limiter
         real(dp) :: Riemann_shock_D_mix_reduction_on, &
            Riemann_shock_D_mix_reduction_full_on
         logical :: use_face_reconstruction
         logical :: include_mlt_in_velocity_time_centering
         logical :: use_hydro_merge_limits_in_mesh_plan
         logical :: split_merge_amr_use_metric_zoning_for_u_flag
         logical :: split_merge_amr_reconstruct_pressure_for_u_flag
         real(dp) :: split_merge_amr_metric_logR_weight, &
            split_merge_amr_metric_logtau_weight, &
            split_merge_amr_metric_min_delta_lnR, &
            split_merge_amr_metric_min_delta_lntau

         real(dp) :: RSP2_Lsurf_factor
         real(dp) :: RSP2_alfap
         real(dp) :: RSP2_alfat
         real(dp) :: RSP2_alfam
         real(dp) :: RSP2_alfar
         real(dp) :: RSP2_alfad

         real(dp) :: RSP2_min_Lt_div_L_for_overshooting_mixing_type
         real(dp) :: RSP2_min_Lc_div_L_for_convective_mixing_type

         real(dp) :: RSP2_w_fix_if_neg, RSP2_source_seed

         logical :: RSP2_use_L_eqn_at_surface, RSP2_assume_HSE, RSP2_report_adjust_w
         logical :: RSP2_use_mass_interp_face_values
         logical :: RSP2_use_Stellingwerf_Lr
         integer :: RSP2_num_outermost_cells_forced_nonturbulent, &
            RSP2_num_innermost_cells_forced_nonturbulent

         logical :: RSP2_remesh_when_load
         integer :: RSP2_nz_div_IBOTOM

         integer :: RSP2_target_steps_per_cycle, RSP2_max_num_periods, RSP2_work_period, &
            RSP2_map_first_period, RSP2_map_last_period, RSP2_map_zone_interval
         real(dp) :: RSP2_min_max_R_for_periods, RSP2_GREKM_avg_abs_frac_new, &
            RSP2_GREKM_avg_abs_limit
         character(256) :: RSP2_work_filename, RSP2_map_columns_filename, &
            RSP2_map_filename, RSP2_map_history_filename
         logical :: RSP2_write_map

         logical :: star_LNA_flag
         integer :: star_LNA_model_number
         logical :: star_LNA_stop_after_run
         real(dp) :: star_LNA_T_inner
         integer :: star_LNA_num_modes
         real(dp) :: star_LNA_min_mode_frequency_uHz
         real(dp) :: star_LNA_min_first_mode_eta
         real(dp) :: star_LNA_min_mode_eta
         real(dp) :: star_LNA_max_abs_mode_eta
         real(dp) :: star_LNA_max_eigenvector_residual
         character(256) :: star_LNA_output_directory
         character(256) :: star_LNA_output_file_prefix
         logical :: star_LNA_write_matrix_summary
         logical :: star_LNA_write_period_growth
         logical :: star_LNA_write_eigenfunctions
         logical :: star_LNA_write_work_integrals
         logical :: star_LNA_set_initial_velocity
         real(dp) :: star_LNA_kick_vsurf_km_per_sec
         integer :: star_LNA_mode_for_period
         integer :: star_LNA_kick_mode_1, star_LNA_kick_mode_2, star_LNA_kick_mode_3
         real(dp) :: star_LNA_kick_fraction_1, star_LNA_kick_fraction_2, star_LNA_kick_fraction_3
         character(32) :: star_LNA_solver
         logical :: star_LNA_include_tdc, star_LNA_include_rsp2
         logical :: star_LNA_perturb_convective_flux
         logical :: star_LNA_perturb_turbulent_pressure
         logical :: star_LNA_perturb_eddy_viscosity
         logical :: star_LNA_perturb_turbulent_energy
         character(32) :: star_LNA_convection_treatment

         logical :: gyre_write_tdc_lna_background



         logical :: use_mass_corrections

         real(dp) :: op_split_burn_min_T_for_variable_T_solver

      ! WD phase separation

         logical :: phase_separation_no_diffusion
