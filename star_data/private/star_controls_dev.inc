         real(dp) :: TDC_Lsurf_factor
         real(dp) :: TDC_alfa
         real(dp) :: TDC_alfap
         real(dp) :: TDC_alfat
         real(dp) :: TDC_alfam
         real(dp) :: TDC_alfar
         
         real(dp) :: TDC_min_Lt_div_L_for_overshooting_mixing_type
         real(dp) :: TDC_min_Lc_div_L_for_convective_mixing_type
         
         logical :: TDC_use_L_eqn_at_surface, TDC_assume_HSE
         logical :: TDC_use_RSP_eqn_for_Y_face
         logical :: TDC_use_mass_interp_face_values
         logical :: TDC_use_Stellingwerf_Lr
         real(dp) :: TDC_w_min_for_damping, TDC_source_seed
         integer :: TDC_num_outermost_cells_forced_nonturbulent, &
            TDC_num_innermost_cells_forced_nonturbulent

         integer :: TDC_target_steps_per_cycle, TDC_max_num_periods, TDC_work_period, &
            TDC_map_first_period, TDC_map_last_period, TDC_map_zone_interval
         real(dp) :: TDC_min_max_R_for_periods, TDC_GREKM_avg_abs_frac_new, &
            TDC_GREKM_avg_abs_limit
         character(256) :: TDC_work_filename, TDC_map_columns_filename, &
            TDC_map_filename, TDC_map_history_filename
         logical :: TDC_write_map

            
            
         
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
         