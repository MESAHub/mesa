! star_data_procedures.inc

         procedure(id_only_function_interface), pointer, nopass :: &
            finish_relax_step => null()
         procedure(id_only_subroutine_interface), pointer, nopass :: &
            finished_relax => null()
           
         procedure(how_many_extra_header_items_interface), pointer, nopass :: &
            how_many_extra_profile_header_items => null()
         procedure(data_for_extra_header_items_interface), pointer, nopass :: &
            data_for_extra_profile_header_items => null()
           
         procedure(how_many_extra_header_items_interface), pointer, nopass :: &
            how_many_extra_history_header_items => null()
         procedure(data_for_extra_header_items_interface), pointer, nopass :: &
            data_for_extra_history_header_items => null()
            
      ! "extras" procedures
         procedure(extras_startup_interface), pointer, nopass :: &
            extras_startup => null()
            
         procedure(extras_check_model_interface), pointer, nopass :: &
            extras_check_model => null()
            
         procedure(extras_start_step_interface), pointer, nopass :: &
            extras_start_step => null()
            
         procedure(extras_finish_step_interface), pointer, nopass :: &
            extras_finish_step => null()
            
         procedure(extras_after_evolve_interface), pointer, nopass :: &
            extras_after_evolve => null()
         
         
      ! "other" procedures

         procedure(other_remove_surface_interface), pointer, nopass :: &
            other_remove_surface => null()         
         procedure(other_adjust_mdot_interface), pointer, nopass :: &
            other_adjust_mdot => null()
         procedure(other_j_for_adjust_J_lost_interface), pointer, nopass :: &
            other_j_for_adjust_J_lost => null()
         procedure(other_after_enter_setmatrix_interface), pointer, nopass :: &
            other_after_enter_setmatrix => null()
         procedure(other_after_struct_burn_mix_interface), pointer, nopass :: &
            other_after_struct_burn_mix => null()
         procedure(other_before_struct_burn_mix_interface), pointer, nopass :: &
            other_before_struct_burn_mix => null()
         procedure(other_wind_interface), pointer, nopass :: &
            other_wind => null()
         procedure(other_accreting_state_interface), pointer, nopass :: &
            other_accreting_state => null()
         procedure (other_D_mix_interface), pointer, nopass :: &
            other_D_mix => null()
         procedure (other_am_mixing_interface), pointer, nopass :: &
            other_am_mixing => null()
         procedure (other_momentum_implicit_interface), pointer, nopass :: &
            other_momentum_implicit => null()
         procedure (other_momentum_interface), pointer, nopass :: &
            other_momentum => null()
         procedure (other_pressure_interface), pointer, nopass :: &
            other_pressure => null()
         procedure (other_energy_interface), pointer, nopass :: &
            other_energy => null()
         procedure (other_diffusion_coefficients_interface), pointer, nopass :: &
            other_diffusion_coefficients => null()
         procedure (other_energy_implicit_interface), pointer, nopass :: &
            other_energy_implicit => null()
         procedure (other_eval_fp_ft_interface), pointer, nopass :: &
            other_eval_fp_ft => null()            
         procedure (other_eval_i_rot_interface), pointer, nopass :: &
            other_eval_i_rot => null()            
         procedure (other_torque_interface), pointer, nopass :: &
            other_torque => null()            
         procedure (other_torque_implicit_interface), pointer, nopass :: &
            other_torque_implicit => null()            
         procedure (other_brunt_interface), pointer, nopass :: &
            other_brunt => null()            
         procedure (other_brunt_smoothing_interface), pointer, nopass :: &
            other_brunt_smoothing => null()            
         procedure (other_build_initial_model_interface), pointer, nopass :: &
            other_build_initial_model => null()            
         procedure (other_cgrav_interface), pointer, nopass :: &
            other_cgrav => null()
         procedure (other_mesh_delta_coeff_factor_interface), pointer, nopass :: &
            other_mesh_delta_coeff_factor => null()
         procedure (other_alpha_mlt_interface), pointer, nopass :: &
            other_alpha_mlt => null()
         procedure (other_opacity_factor_interface), pointer, nopass :: &
            other_opacity_factor => null()
         procedure (other_diffusion_factor_interface), pointer, nopass :: &
            other_diffusion_factor => null()
         procedure (other_eps_grav_interface), pointer, nopass :: &
            other_eps_grav => null()
         procedure (other_rsp_linear_analysis_interface), pointer, nopass :: &
            other_rsp_linear_analysis => null()
         procedure (other_rsp_build_model_interface), pointer, nopass :: &
            other_rsp_build_model => null()
         procedure (other_gradr_factor_interface), pointer, nopass :: &
            other_gradr_factor => null()
         procedure (other_overshooting_scheme_interface), pointer, nopass :: &
            other_overshooting_scheme => null()
         
         procedure (other_mlt_interface), pointer, nopass :: &
            other_mlt => null()
         procedure (other_solver_monitor_interface), pointer, nopass :: &
            other_solver_monitor => null()
         procedure (other_neu_interface), pointer, nopass :: &
            other_neu => null()
         procedure (other_net_get_interface), pointer, nopass :: &
            other_net_get => null()
         procedure (other_adjust_mlt_gradT_fraction_interface), pointer, nopass :: &
            other_adjust_mlt_gradT_fraction => null()
         procedure (other_after_set_mixing_info_interface), pointer, nopass :: &
            other_after_set_mixing_info => null()
         procedure (other_diffusion_interface), pointer, nopass :: &
            other_diffusion => null()

         procedure (other_surface_PT_interface), pointer, nopass :: &
            other_surface_PT => null()

         procedure (other_eosDT_get_interface), pointer, nopass :: &
            other_eosDT_get => null()
         procedure (other_eosDT_get_T_interface), pointer, nopass :: &
            other_eosDT_get_T => null()
         procedure (other_eosDT_get_Rho_interface), pointer, nopass :: &
            other_eosDT_get_Rho => null()

         procedure (other_eosPT_get_interface), pointer, nopass :: &
            other_eosPT_get => null()
         procedure (other_eosPT_get_T_interface), pointer, nopass :: &
            other_eosPT_get_T => null()
         procedure (other_eosPT_get_Pgas_interface), pointer, nopass :: &
            other_eosPT_get_Pgas => null()
         procedure (other_eosPT_get_Pgas_for_Rho_interface), pointer, nopass :: &
            other_eosPT_get_Pgas_for_Rho => null()

         procedure (other_eosDE_get_interface), pointer, nopass :: &
            other_eosDE_get => null()

         procedure (other_kap_get_interface), pointer, nopass :: &
            other_kap_get => null()
         procedure (other_kap_get_op_mono_interface), pointer, nopass :: &
            other_kap_get_op_mono => null()

         procedure (other_pgstar_plots_info_interface), pointer, nopass :: &
            other_pgstar_plots_info => null()
         procedure (how_many_other_mesh_fcns_interface), pointer, nopass :: &
            how_many_other_mesh_fcns => null()
         procedure (other_mesh_fcn_data_interface), pointer, nopass :: &
            other_mesh_fcn_data => null()

         procedure (other_write_pulsation_info_interface), pointer, nopass :: &
            other_write_pulsation_info => null()
         procedure (other_get_pulsation_info_interface), pointer, nopass :: &
            other_get_pulsation_info => null()

         procedure (other_photo_write_interface), pointer, nopass :: &
            other_photo_write => null()
         procedure (other_photo_read_interface), pointer, nopass :: &
            other_photo_read => null()

         procedure (other_new_generation_interface), pointer, nopass :: &
            other_new_generation => null()
         procedure (other_set_current_to_old_interface), pointer, nopass :: &
            other_set_current_to_old => null()

         procedure (other_astero_freq_corr_interface), pointer, nopass :: &
            other_astero_freq_corr => null()

         procedure (other_timestep_limit_interface), pointer, nopass :: &
            other_timestep_limit => null()

         procedure (how_many_extra_history_columns_interface), pointer, nopass :: &
            how_many_extra_history_columns => null()
         procedure (data_for_extra_history_columns_interface), pointer, nopass :: &
            data_for_extra_history_columns => null()

         procedure (how_many_extra_profile_columns_interface), pointer, nopass :: &
            how_many_extra_profile_columns => null()
         procedure (data_for_extra_profile_columns_interface), pointer, nopass :: &
            data_for_extra_profile_columns => null()

         procedure (how_many_binary_history_columns_interface), pointer, nopass :: &
            how_many_binary_history_columns => null()
         procedure (data_for_binary_history_columns_interface), pointer, nopass :: &
            data_for_binary_history_columns => null()

         procedure (how_many_extra_binary_history_columns_interface), pointer, nopass :: &
            how_many_extra_binary_history_columns => null()
         procedure (data_for_extra_binary_history_columns_interface), pointer, nopass :: &
            data_for_extra_binary_history_columns => null()

         procedure(other_export_pulse_data_interface), pointer, nopass :: &
              other_export_pulse_data => null()
         procedure(other_get_pulse_data_interface), pointer, nopass :: &
              other_get_pulse_data => null()
         procedure (other_edit_pulse_data_interface), pointer, nopass :: &
              other_edit_pulse_data => null()

         procedure (binary_other_torque_implicit_interface), pointer, nopass :: &
            binary_other_torque_implicit => null()
         procedure (binary_other_torque_interface), pointer, nopass :: &
            binary_other_torque => null()

         procedure(pgstar_decorator_interface), pointer, nopass :: &
            Abundance_pgstar_decorator => null(), &
            Color_Magnitude1_pgstar_decorator => null(), &
            Color_Magnitude2_pgstar_decorator => null(), &
            Color_Magnitude3_pgstar_decorator => null(), &
            Color_Magnitude4_pgstar_decorator => null(), &
            Color_Magnitude5_pgstar_decorator => null(), &
            Color_Magnitude6_pgstar_decorator => null(), &
            Color_Magnitude7_pgstar_decorator => null(), &
            Color_Magnitude8_pgstar_decorator => null(), &
            Color_Magnitude9_pgstar_decorator => null(), &
            dPg_dnu_pgstar_decorator => null(), &
            Dynamo_pgstar_decorator => null(), &
            History_Panels1_pgstar_decorator => null(), &
            History_Panels2_pgstar_decorator => null(), &
            History_Panels3_pgstar_decorator => null(), &
            History_Panels4_pgstar_decorator => null(), &
            History_Panels5_pgstar_decorator => null(), &
            History_Panels6_pgstar_decorator => null(), &
            History_Panels7_pgstar_decorator => null(), &
            History_Panels8_pgstar_decorator => null(), &
            History_Panels9_pgstar_decorator => null(), &
            History_Track1_pgstar_decorator => null(), &
            History_Track2_pgstar_decorator => null(), &
            History_Track3_pgstar_decorator => null(), &
            History_Track4_pgstar_decorator => null(), &
            History_Track5_pgstar_decorator => null(), &
            History_Track6_pgstar_decorator => null(), &
            History_Track7_pgstar_decorator => null(), &
            History_Track8_pgstar_decorator => null(), &
            History_Track9_pgstar_decorator => null(), &
            HR_pgstar_decorator => null(), &
            Kipp_pgstar_decorator => null(), &
            logg_logt_pgstar_decorator => null(), &
            logg_teff_pgstar_decorator => null(), &
            logl_r_pgstar_decorator => null(), &
            logl_teff_pgstar_decorator => null(), &
            logl_v_pgstar_decorator => null(), &
            l_r_pgstar_decorator => null(), &
            l_teff_pgstar_decorator => null(), &
            l_v_pgstar_decorator => null(), &
            mixing_pgstar_decorator => null(), &
            mode_prop_pgstar_decorator => null(), &
            network_pgstar_decorator => null(), &
            power_pgstar_decorator => null(), &
            production_pgstar_decorator => null(), &
            Profile_Panels1_pgstar_decorator => null(), &
            Profile_Panels2_pgstar_decorator => null(), &
            Profile_Panels3_pgstar_decorator => null(), &
            Profile_Panels4_pgstar_decorator => null(), &
            Profile_Panels5_pgstar_decorator => null(), &
            Profile_Panels6_pgstar_decorator => null(), &
            Profile_Panels7_pgstar_decorator => null(), &
            Profile_Panels8_pgstar_decorator => null(), &
            Profile_Panels9_pgstar_decorator => null(), &
            R_L_pgstar_decorator => null(), &
            R_Teff_pgstar_decorator => null(), &
            rti_pgstar_decorator => null(), &
            summary_burn_pgstar_decorator => null(), &
            summary_profile_pgstar_decorator => null(), &
            summary_history_pgstar_decorator => null(), &
            trho_pgstar_decorator => null(), &
            tmaxrho_pgstar_decorator => null(), &
            trho_profile_pgstar_decorator => null()
