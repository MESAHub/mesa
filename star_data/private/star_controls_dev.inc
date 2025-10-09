         
         logical :: compare_TDC_to_MLT 

         real(dp) :: alpha_TDC_DAMPM
         real(dp) :: alpha_TDC_C
         real(dp) :: alpha_TDC_S
         logical :: remesh_for_TDC_pulsations_log_core_zoning
         logical :: TDC_use_density_form_for_eddy_viscosity
         logical :: use_RSP_L_eqn_outer_BC
         real(dp) :: TDC_num_innermost_cells_forced_nonturbulent
         logical :: include_mlt_Pturb_in_thermodynamic_gradients
         logical :: make_mlt_hydrodynamic, include_mlt_corr_to_TDC
         logical :: TDC_include_eturb_in_energy_equation, use_rsp_form_of_scale_height
         logical :: include_mlt_in_velocity_time_centering
         logical :: TDC_hydro_use_mass_interp_face_values
         integer :: TDC_hydro_nz, TDC_hydro_nz_outer
         real(dp) :: TDC_hydro_T_anchor, TDC_hydro_dq_1_factor

         real(dp) :: RSP2_Lsurf_factor
         real(dp) :: RSP2_alfap
         real(dp) :: RSP2_alfat
         real(dp) :: RSP2_alfam
         real(dp) :: RSP2_alfar
         real(dp) :: RSP2_alfad

         real(dp) :: RSP2_min_Lt_div_L_for_overshooting_mixing_type
         real(dp) :: RSP2_min_Lc_div_L_for_convective_mixing_type

         real(dp) :: RSP2_w_fix_if_neg, RSP2_w_min_for_damping, RSP2_source_seed

         logical :: RSP2_use_L_eqn_at_surface, RSP2_assume_HSE, RSP2_report_adjust_w
         logical :: RSP2_use_RSP_eqn_for_Y_face
         logical :: RSP2_use_mass_interp_face_values
         logical :: RSP2_use_Stellingwerf_Lr
         integer :: RSP2_num_outermost_cells_forced_nonturbulent, &
            RSP2_num_innermost_cells_forced_nonturbulent

         logical :: RSP2_remesh_when_load
         integer :: RSP2_nz, RSP2_nz_outer, RSP2_nz_div_IBOTOM
         real(dp) :: RSP2_T_anchor, RSP2_dq_1_factor

         integer :: RSP2_target_steps_per_cycle, RSP2_max_num_periods, RSP2_work_period, &
            RSP2_map_first_period, RSP2_map_last_period, RSP2_map_zone_interval
         real(dp) :: RSP2_min_max_R_for_periods, RSP2_GREKM_avg_abs_frac_new, &
            RSP2_GREKM_avg_abs_limit
         character(256) :: RSP2_work_filename, RSP2_map_columns_filename, &
            RSP2_map_filename, RSP2_map_history_filename
         logical :: RSP2_write_map



         logical :: use_mass_corrections

         real(dp) :: op_split_burn_min_T_for_variable_T_solver

      ! WD phase separation

         logical :: phase_separation_no_diffusion
