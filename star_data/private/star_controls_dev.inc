
         real(dp) :: max_dt_div_tau_conv_for_TDC, max_dt_years_for_TDC
         real(dp) :: alpha_TDC_DAMP   ! ~ RSP2_alfad
         real(dp) :: alpha_TDC_DAMPR  ! ~ RSP2_alfar
         real(dp) :: alpha_TDC_eddy_viscosity   ! RSP2_alfam
         real(dp) :: alpha_TDC_turbulent_pressure  ! RSP2_alfap
         
         real(dp) :: max_X_for_TDC
         real(dp) :: max_X_for_gradT_eqn
         logical :: compare_TDC_to_MLT





         real(dp) :: RSP2_Lsurf_factor
         real(dp) :: RSP2_alfap
         real(dp) :: RSP2_alfat
         real(dp) :: RSP2_alfam
         real(dp) :: RSP2_alfar
         real(dp) :: RSP2_alfad
         
         real(dp) :: RSP2_min_Lt_div_L_for_overshooting_mixing_type
         real(dp) :: RSP2_min_Lc_div_L_for_convective_mixing_type
         
         real(dp) :: RSP2_min_dt_div_tau_conv_switch_to_MLT
         real(dp) :: RSP2_min_dt_years_switch_to_MLT
         
         real(dp) :: RSP2_w_min_for_damping, RSP2_source_seed
         
         logical :: RSP2_use_L_eqn_at_surface, RSP2_assume_HSE
         logical :: RSP2_use_RSP_eqn_for_Y_face, include_w_in_correction_limits
         logical :: RSP2_use_mass_interp_face_values
         logical :: RSP2_use_Stellingwerf_Lr
         integer :: RSP2_num_outermost_cells_forced_nonturbulent, &
            RSP2_num_innermost_cells_forced_nonturbulent

         integer :: RSP2_target_steps_per_cycle, RSP2_max_num_periods, RSP2_work_period, &
            RSP2_map_first_period, RSP2_map_last_period, RSP2_map_zone_interval
         real(dp) :: RSP2_min_max_R_for_periods, RSP2_GREKM_avg_abs_frac_new, &
            RSP2_GREKM_avg_abs_limit
         character(256) :: RSP2_work_filename, RSP2_map_columns_filename, &
            RSP2_map_filename, RSP2_map_history_filename
         logical :: RSP2_write_map

            
            
         
         logical :: use_superad_reduction
         real(dp) :: &
            superad_reduction_Gamma_limit, &
            superad_reduction_Gamma_limit_scale, &
            superad_reduction_Gamma_inv_scale, &
            superad_reduction_diff_grads_limit, &
            superad_reduction_limit
            
            

         real(dp) :: &
            conv_vel_mix_factor, &
            conv_vel_drag_factor, &
            conv_vel_buoyancy_factor, &
            conv_vel_D, &
            conv_vel_siglimit, &
            conv_vel_v0
            
         real(dp) :: min_q_for_normal_mlt_gradT_full_off, max_q_for_normal_mlt_gradT_full_on
         
         logical :: conv_vel_ignore_thermohaline, conv_vel_ignore_semiconvection, &
                    conv_vel_fully_lagrangian, conv_vel_include_homologous_term, &
                    conv_vel_use_mlt_vc_start

                    !convergence_separate_equ_conv_vel_residuals, &
                    !convergence_separate_equ_conv_vel_corrections, &

         !real(dp) :: &
         !   conv_vel_tol_max_residual, &
         !   conv_vel_tol_max_correction, &
         !   conv_vel_min_correction_to_ignore_residual, &
         !   conv_vel_scale_max_correction, &
         !   conv_vel_corr_coeff_limit
         !integer :: conv_vel_solver_itermin_until_reduce_min_corr_coeff
         