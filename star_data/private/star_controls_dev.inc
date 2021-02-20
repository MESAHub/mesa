      
         real(dp) :: TDC_Lsurf_factor
         real(dp) :: TDC_alfa
         real(dp) :: TDC_alfap
         real(dp) :: TDC_alfat
         real(dp) :: TDC_alfam
         real(dp) :: TDC_alfar
         
         logical :: TDC_use_Stellingwerf_Lr
         real(dp) :: TDC_w_min_for_damping
         integer :: TDC_num_outermost_cells_forced_nonturbulent, &
            TDC_num_innermost_cells_forced_nonturbulent
         
         
         real(dp) :: accretion_entropy_at_surface
         
         real(dp) :: D_mix_zero_region_bottom_q, D_mix_zero_region_top_q, &
            dq_D_mix_zero_at_H_He_crossover, dq_D_mix_zero_at_H_C_crossover
                  

         logical :: use_DGESVX_in_bcyclic, use_equilibration_in_DGESVX
         logical :: report_min_rcond_from_DGESXV
         
         

         logical :: mlt_make_surface_no_mixing
         
         logical :: do_normalize_dqs_as_part_of_set_qs
            
            
         
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
         