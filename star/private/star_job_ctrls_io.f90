! ***********************************************************************
!
!   Copyright (C) 2013  The MESA Team
!
!   MESA is free software; you can use it and/or modify
!   it under the combined terms and restrictions of the MESA MANIFESTO
!   and the GNU General Library Public License as published
!   by the Free Software Foundation; either version 2 of the License,
!   or (at your option) any later version.
!
!   You should have received a copy of the MESA MANIFESTO along with
!   this software; if not, it is available at the mesa website:
!   http://mesa.sourceforge.net/
!
!   MESA is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!   See the GNU Library General Public License for more details.
!
!   You should have received a copy of the GNU Library General Public License
!   along with this software; if not, write to the Free Software
!   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
!
! ***********************************************************************

module star_job_ctrls_io
   
   use const_def, only : dp
   use star_private_def
   use rates_def, only : maxlen_reaction_Name
   
   implicit none
   
   include 'star_job_controls.inc'
   include 'star_job_controls_dev.inc'
   
   namelist /star_job/ &
      mesa_dir, &
      eosDT_cache_dir, &
      kap_cache_dir, &
      rates_cache_dir, &
      pause_before_terminate, &
      profile_columns_file, &
      history_columns_file, &
      show_log_description_at_start, &
      list_net_reactions, &
      show_net_reactions_info, &
      show_net_species_info, &
      first_model_for_timing, &
      steps_before_start_timing, &
      show_eqns_and_vars_names, &
      pgstar_flag, &
      disable_pgstar_during_relax_flag, &
      clear_initial_pgstar_history, &
      clear_pgstar_history, &
      save_pgstar_files_when_terminate, &
      save_photo_when_terminate, &
      load_saved_photo, &
      saved_photo_name, &
      save_photo_filename, &
      save_photo_number, &
      write_profile_when_terminate, &
      show_retry_counts_when_terminate, &
      show_timestep_limit_counts_when_terminate, &
      filename_for_profile_when_terminate, &
      create_pre_main_sequence_model, &
      pre_ms_relax_to_start_radiative_core, &
      pre_ms_relax_num_steps, &
      pre_ms_min_steps_before_check_radiative_core, &
      pre_ms_check_radiative_core_start, &
      pre_ms_check_radiative_core_stop, &
      pre_ms_check_radiative_core_Lnuc_div_L_limit, &
      pre_ms_check_radiative_core_min_mass, &
      pre_ms_T_c, &
      pre_ms_guess_rho_c, &
      pre_ms_d_log10_P, &
      pre_ms_logT_surf_limit, &
      pre_ms_logP_surf_limit, &
      create_initial_model, &
      initial_model_relax_num_steps, &
      radius_in_cm_for_create_initial_model, &
      mass_in_gm_for_create_initial_model, &
      center_logP_1st_try_for_create_initial_model, &
      entropy_1st_try_for_create_initial_model, &
      abs_e01_tolerance_for_create_initial_model, &
      abs_e02_tolerance_for_create_initial_model, &
      max_tries_for_create_initial_model, &
      initial_model_eps, &
      save_star_job_namelist, &
      star_job_namelist_name, &
      echo_at_start, &
      echo_at_end, &
      load_saved_model, &
      load_model_filename, &
      create_merger_model, &
      saved_model_for_merger_1, &
      saved_model_for_merger_2, &
      set_max_dt_to_frac_lifetime, &
      max_frac_of_lifetime_per_step, &
      astero_just_call_my_extras_check_model, &
      relax_mass, &
      relax_initial_mass, &
      new_mass, &
      lg_max_abs_mdot, &
      relax_mass_to_remove_H_env, &
      relax_initial_mass_to_remove_H_env, &
      extra_mass_retained_by_remove_H_env, &
      relax_mass_scale, &
      relax_initial_mass_scale, &
      dlgm_per_step, &
      
      relax_M_center_dt, &
      change_mass_years_for_dt, &
      relax_M_center, &
      relax_initial_M_center, &
      relax_core, &
      relax_initial_core, &
      new_core_mass, &
      dlg_core_mass_per_step, &
      relax_core_years_for_dt, &
      core_avg_rho, &
      core_avg_eps, &
      
      relax_R_center, &
      relax_initial_R_center, &
      new_R_center, &
      dlgR_per_step, &
      relax_R_center_dt, &
      
      set_v_center, &
      set_initial_v_center, &
      relax_v_center, &
      relax_initial_v_center, &
      new_v_center, &
      dv_per_step, &
      relax_v_center_dt, &
      
      zero_alpha_RTI, &
      zero_initial_alpha_RTI, &
      
      set_L_center, &
      set_initial_L_center, &
      relax_L_center, &
      relax_initial_L_center, &
      new_L_center, &
      dlgL_per_step, &
      relax_L_center_dt, &
      
      remove_center_at_cell_k, &
      remove_center_by_temperature, &
      remove_center_by_mass_fraction_q, &
      remove_center_by_delta_mass_gm, &
      remove_center_by_delta_mass_Msun, &
      remove_center_by_mass_gm, &
      remove_center_by_mass_Msun, &
      remove_center_by_radius_cm, &
      remove_center_by_radius_Rsun, &
      remove_center_by_he4, &
      remove_center_by_c12_o16, &
      remove_center_by_si28, &
      remove_center_to_reduce_co56_ni56, &
      remove_center_by_ye, &
      remove_center_by_entropy, &
      remove_center_by_infall_kms, &
      remove_center_at_inner_max_abs_v, &
      remove_fe_core, &
      remove_fallback_at_each_step, &
      fallback_check_total_energy, &
      remove_fallback_speed_limit, &
      remove_center_set_zero_v_center, &
      retain_fallback_at_each_step, &
      limit_center_logP_at_each_step, &
      remove_center_adjust_L_center, &
      remove_center_logRho_limit, &
      
      remove_initial_center_at_cell_k, &
      remove_initial_center_by_temperature, &
      remove_initial_center_by_mass_fraction_q, &
      remove_initial_center_by_delta_mass_gm, &
      remove_initial_center_by_delta_mass_Msun, &
      remove_initial_center_by_mass_gm, &
      remove_initial_center_by_mass_Msun, &
      remove_initial_center_by_radius_cm, &
      remove_initial_center_by_radius_Rsun, &
      remove_initial_center_by_he4, &
      remove_initial_center_by_c12_o16, &
      remove_initial_center_by_si28, &
      remove_initial_center_to_reduce_co56_ni56, &
      remove_initial_center_by_ye, &
      remove_initial_center_by_entropy, &
      remove_initial_center_by_infall_kms, &
      remove_initial_center_at_inner_max_abs_v, &
      remove_initial_fe_core, &
      
      zero_initial_inner_v_by_mass_Msun, &
      zero_inner_v_by_mass_Msun, &
      
      remove_surface_at_cell_k, &
      remove_surface_at_he_core_boundary, &
      remove_surface_by_optical_depth, &
      remove_surface_by_density, &
      remove_surface_by_pressure, &
      remove_surface_by_mass_fraction_q, &
      remove_surface_by_mass_gm, &
      remove_surface_by_radius_cm, &
      remove_surface_by_mass_Msun, &
      remove_surface_by_radius_Rsun, &
      remove_surface_by_v_surf_km_s, &
      remove_surface_by_v_surf_div_cs, &
      remove_surface_by_v_surf_div_v_escape, &
      min_q_for_remove_surface_by_v_surf_div_v_escape, &
      max_q_for_remove_surface_by_v_surf_div_v_escape, &
      
      remove_surface_do_jrot, &
      remove_surface_do_entropy, &
      remove_surface_turn_off_energy_sources_and_sinks, &
      remove_surface_by_relax_to_star_cut, &
      
      remove_initial_surface_at_cell_k, &
      remove_initial_surface_at_he_core_boundary, &
      remove_initial_surface_by_optical_depth, &
      remove_initial_surface_by_density, &
      remove_initial_surface_by_pressure, &
      remove_initial_surface_by_mass_fraction_q, &
      remove_initial_surface_by_mass_gm, &
      remove_initial_surface_by_radius_cm, &
      remove_initial_surface_by_mass_Msun, &
      remove_initial_surface_by_radius_Rsun, &
      remove_initial_surface_by_v_surf_km_s, &
      remove_initial_surface_by_v_surf_div_cs, &
      remove_initial_surface_by_v_surf_div_v_escape, &
      
      report_mass_not_fe56, &
      relax_dxdt_nuc_factor, &
      relax_initial_dxdt_nuc_factor, &
      new_dxdt_nuc_factor, &
      dxdt_nuc_factor_multiplier, &
      relax_eps_nuc_factor, &
      relax_initial_eps_nuc_factor, &
      new_eps_nuc_factor, &
      eps_nuc_factor_multiplier, &
      relax_opacity_max, &
      relax_initial_opacity_max, &
      new_opacity_max, &
      opacity_max_multiplier, &
      relax_max_surf_dq, &
      relax_initial_max_surf_dq, &
      new_max_surf_dq, &
      max_surf_dq_multiplier, &
      
      relax_tau_factor, &
      relax_initial_tau_factor, &
      set_tau_factor, &
      set_initial_tau_factor, &
      relax_to_this_tau_factor, &
      set_to_this_tau_factor, &
      dlogtau_factor, &
      set_tau_factor_after_core_He_burn, &
      set_tau_factor_after_core_C_burn, &
      relax_tau_factor_after_core_He_burn, &
      relax_tau_factor_after_core_C_burn, &
      
      relax_to_this_opacity_factor, &
      d_opacity_factor, &
      relax_opacity_factor, &
      relax_initial_opacity_factor, &
      
      adjust_tau_factor_to_surf_density, &
      base_for_adjust_tau_factor_to_surf_density, &
      
      relax_Tsurf_factor, &
      relax_initial_Tsurf_factor, &
      set_Tsurf_factor, &
      set_initial_Tsurf_factor, &
      relax_to_this_Tsurf_factor, &
      set_to_this_Tsurf_factor, &
      dlogTsurf_factor, &
      
      relax_irradiation, &
      relax_initial_irradiation, &
      set_irradiation, &
      set_initial_irradiation, &
      relax_irradiation_min_steps, &
      relax_to_this_irrad_flux, &
      set_to_this_irrad_flux, &
      irrad_col_depth, &
      relax_irradiation_max_yrs_dt, &
      relax_mass_change, &
      relax_initial_mass_change, &
      relax_mass_change_min_steps, &
      relax_mass_change_max_yrs_dt, &
      relax_mass_change_init_mdot, &
      relax_mass_change_final_mdot, &
      
      change_RTI_flag, &
      change_initial_RTI_flag, &
      new_RTI_flag, &
      
      change_RSP_flag, &
      change_initial_RSP_flag, &
      new_RSP_flag, &
      
      change_RSP2_flag, &
      change_initial_RSP2_flag, &
      change_RSP2_flag_at_model_number, &
      new_RSP2_flag, &
      create_RSP2_model, &
      
      change_w_div_wc_flag, &
      change_initial_w_div_wc_flag, &
      new_w_div_wc_flag, &
      
      change_j_rot_flag, &
      change_initial_j_rot_flag, &
      new_j_rot_flag, &
      
      create_RSP_model, &
      
      change_v_flag, &
      change_initial_v_flag, &
      new_v_flag, &
      
      change_D_omega_flag, &
      change_initial_D_omega_flag, &
      new_D_omega_flag, &
      
      change_am_nu_rot_flag, &
      change_initial_am_nu_rot_flag, &
      new_am_nu_rot_flag, &
      
      use_D_omega_for_am_nu_rot, &
      
      change_u_flag, &
      change_initial_u_flag, &
      new_u_flag, &
      
      change_reconstruction_flag, &
      change_initial_reconstruction_flag, &
      new_reconstruction_flag, &
      
      center_ye_limit_for_v_flag, &
      change_rotation_flag, &
      change_initial_rotation_flag, &
      new_rotation_flag, &
      use_w_div_wc_flag_with_rotation, &
      use_j_rot_flag_with_rotation, &
      
      
      set_omega, &
      set_initial_omega, &
      set_omega_step_limit, &
      set_near_zams_omega_steps, &
      new_omega, &
      set_omega_div_omega_crit, &
      set_initial_omega_div_omega_crit, &
      set_omega_div_omega_crit_step_limit, &
      set_near_zams_omega_div_omega_crit_steps, &
      new_omega_div_omega_crit, &
      set_surface_rotation_v, &
      set_initial_surface_rotation_v, &
      set_surf_rotation_v_step_limit, &
      set_near_zams_surface_rotation_v_steps, &
      new_surface_rotation_v, &
      relax_omega, &
      relax_initial_omega, &
      near_zams_relax_omega, &
      relax_omega_div_omega_crit, &
      relax_initial_omega_div_omega_crit, &
      near_zams_relax_omega_div_omega_crit, &
      relax_surface_rotation_v, &
      relax_initial_surface_rotation_v, &
      near_zams_relax_initial_surface_rotation_v, &
      num_steps_to_relax_rotation, &
      relax_omega_max_yrs_dt, &
      set_uniform_initial_composition, &
      initial_h1, &
      initial_h2, &
      initial_he3, &
      initial_he4, &
      initial_zfracs, &
      dump_missing_metals_into_heaviest, &
      relax_initial_composition, &
      relax_initial_to_xaccrete, &
      relax_composition_filename, &
      num_steps_to_relax_composition, &
      timescale_for_relax_composition, &
      relax_initial_angular_momentum, &
      max_steps_to_relax_angular_momentum, &
      timescale_for_relax_angular_momentum, &
      max_dt_for_relax_angular_momentum, &
      relax_angular_momentum_constant_omega_center, &
      num_timescales_for_relax_angular_momentum, &
      relax_angular_momentum_filename, &
      relax_initial_entropy, &
      max_steps_to_relax_entropy, &
      timescale_for_relax_entropy, &
      max_dt_for_relax_entropy, &
      num_timescales_for_relax_entropy, &
      relax_entropy_filename, &
      get_entropy_for_relax_from_eos, &
      report_cell_for_xm, &
      set_to_xa_for_accretion, &
      set_initial_to_xa_for_accretion, &
      set_nzlo, &
      set_nzhi, &
      change_Y, &
      change_initial_Y, &
      relax_Y, &
      relax_initial_Y, &
      relax_Y_minq, &
      relax_Y_maxq, &
      new_Y, &
      change_Z, &
      change_initial_Z, &
      relax_Z, &
      relax_initial_Z, &
      relax_Z_minq, &
      relax_Z_maxq, &
      new_Z, &
      steps_to_take_before_terminate, &
      stop_if_this_file_exists, &
      set_initial_age, &
      initial_age, &
      set_initial_model_number, &
      initial_model_number, &
      set_initial_number_retries, &
      initial_number_retries, &
      set_initial_dt, &
      limit_initial_dt, &
      years_for_initial_dt, &
      seconds_for_initial_dt, &
      
      set_initial_cumulative_energy_error, &
      set_cumulative_energy_error, &
      set_cumulative_energy_error_at_step, &
      set_cumulative_energy_error_each_step_if_age_less_than, &
      set_cumulative_energy_error_each_relax, &
      new_cumulative_energy_error, &
      
      change_net, &
      change_initial_net, &
      new_net_name, &
      change_small_net, &
      change_initial_small_net, &
      new_small_net_name, &
      
      h_he_net, &
      co_net, &
      adv_net, &
      adjust_abundances_for_new_isos, &
      set_uniform_xa_from_file, &
      set_uniform_initial_xa_from_file, &
      file_for_uniform_xa, &
      
      mix_section, mix_initial_section, &
      mix_section_nzlo, mix_section_nzhi, &
      
      T9_weaklib_full_off, &
      T9_weaklib_full_on, &
      weaklib_blend_hi_Z, &
      T9_weaklib_full_off_hi_Z, &
      T9_weaklib_full_on_hi_Z, &
      
      use_suzuki_weak_rates, &
      use_3a_fl87, &
      
      use_special_weak_rates, &
      special_weak_states_file, &
      special_weak_transitions_file, &
      ion_coulomb_corrections, &
      electron_coulomb_corrections, &
      
      mix_envelope_down_to_T, &
      mix_initial_envelope_down_to_T, &
      auto_extend_net, &
      
      enable_adaptive_network, &
      min_x_for_keep, &
      min_x_for_n, &
      min_x_for_add, &
      max_Z_for_add, &
      max_N_for_add, &
      max_A_for_add, &
      
      save_model_number, &
      save_model_filename, &
      save_model_when_terminate, &
      required_termination_code_string, &
      profile_starting_model, &
      profile_model_number, &
      report_retries, &
      
      net_reaction_filename, &
      jina_reaclib_filename, &
      jina_reaclib_min_T9, &
      rate_tables_dir, &
      rate_cache_suffix, &
      read_extra_star_job_inlist1, &
      extra_star_job_inlist1_name, &
      read_extra_star_job_inlist2, &
      extra_star_job_inlist2_name, &
      read_extra_star_job_inlist3, &
      extra_star_job_inlist3_name, &
      read_extra_star_job_inlist4, &
      extra_star_job_inlist4_name, &
      read_extra_star_job_inlist5, &
      extra_star_job_inlist5_name, &
      set_abundance_nzlo, &
      set_abundance_nzhi, &
      set_abundance, &
      set_initial_abundance, &
      chem_name, &
      new_frac, &
      set_abundance_nzlo, set_abundance_nzhi, &
      replace_element, &
      replace_initial_element, &
      chem_name1, &
      chem_name2, &
      replace_element_nzlo, replace_element_nzhi, &
      do_special_test, &
      
      save_pulse_data_for_model_number, &
      save_pulse_data_when_terminate, &
      save_pulse_data_filename, &
      
      chem_isotopes_filename, &
      extras_lipar, &
      extras_lrpar, &
      extras_lcpar, &
      extras_llpar, &
      extras_ipar, &
      extras_rpar, &
      extras_cpar, &
      extras_lpar, &
      num_special_rate_factors, &
      special_rate_factor, &
      filename_of_special_rate, &
      reaction_for_special_factor, &
      color_num_files, &
      color_file_names, &
      color_num_colors, &
      warn_run_star_extras, &
      
      report_garbage_collection, &
      num_steps_for_garbage_collection

contains
   
   
   subroutine do_read_star_job(s, filename, ierr)
      use star_private_def
      use utils_lib
      type (star_info), pointer :: s
      character(*), intent(in) :: filename
      integer, intent(out) :: ierr
      character (len = strlen) :: star_job_namelist_name
      star_job_namelist_name = ''
      ierr = 0
      call set_default_star_job_controls
      call read_star_job_file(s, filename, 1, ierr)
      call check_star_job_controls(s, ierr)
   end subroutine do_read_star_job
   
   
   recursive subroutine read_star_job_file(s, filename, level, ierr)
      use star_private_def
      use utils_lib
      character(*), intent(in) :: filename
      type (star_info), pointer :: s
      integer, intent(in) :: level
      integer, intent(out) :: ierr
      logical :: read_extra1, read_extra2, read_extra3, read_extra4, read_extra5
      character (len = strlen) :: message, extra1, extra2, extra3, extra4, extra5
      integer :: unit
      
      ierr = 0
      
      if (level >= 10) then
         write(*, *) 'ERROR: too many levels of nested extra star_job inlist files'
         ierr = -1
         return
      end if
      
      if (len_trim(filename) > 0) then
         open(newunit = unit, file = trim(filename), action = 'read', delim = 'quote', status = 'old', iostat = ierr)
         if (ierr /= 0) then
            write(*, *) 'Failed to open control namelist file "' // trim(filename) // '"'
            return
         end if
         read(unit, nml = star_job, iostat = ierr)
         close(unit)
         if (ierr /= 0) then
            write(*, *)
            write(*, *)
            write(*, *)
            write(*, *)
            write(*, '(a)') &
               'Failed while trying to read control namelist file: ' // trim(filename)
            write(*, '(a)') &
               'Perhaps the following runtime error message will help you find the problem.'
            write(*, *)
            open(newunit = unit, file = trim(filename), action = 'read', delim = 'quote', status = 'old', iostat = ierr)
            read(unit, nml = star_job)
            close(unit)
            return
         end if
      end if
      
      call store_star_job_controls(s, ierr)
      
      ! recursive calls to read other inlists
      
      read_extra1 = read_extra_star_job_inlist1
      read_extra_star_job_inlist1 = .false.
      extra1 = extra_star_job_inlist1_name
      extra_star_job_inlist1_name = 'undefined'
      
      read_extra2 = read_extra_star_job_inlist2
      read_extra_star_job_inlist2 = .false.
      extra2 = extra_star_job_inlist2_name
      extra_star_job_inlist2_name = 'undefined'
      
      read_extra3 = read_extra_star_job_inlist3
      read_extra_star_job_inlist3 = .false.
      extra3 = extra_star_job_inlist3_name
      extra_star_job_inlist3_name = 'undefined'
      
      read_extra4 = read_extra_star_job_inlist4
      read_extra_star_job_inlist4 = .false.
      extra4 = extra_star_job_inlist4_name
      extra_star_job_inlist4_name = 'undefined'
      
      read_extra5 = read_extra_star_job_inlist5
      read_extra_star_job_inlist5 = .false.
      extra5 = extra_star_job_inlist5_name
      extra_star_job_inlist5_name = 'undefined'
      
      if (read_extra1) then
         call read_star_job_file(s, extra1, level + 1, ierr)
         if (ierr /= 0) return
      end if
      
      if (read_extra2) then
         call read_star_job_file(s, extra2, level + 1, ierr)
         if (ierr /= 0) return
      end if
      
      if (read_extra3) then
         call read_star_job_file(s, extra3, level + 1, ierr)
         if (ierr /= 0) return
      end if
      
      if (read_extra4) then
         call read_star_job_file(s, extra4, level + 1, ierr)
         if (ierr /= 0) return
      end if
      
      if (read_extra5) then
         call read_star_job_file(s, extra5, level + 1, ierr)
         if (ierr /= 0) return
      end if
   
   end subroutine read_star_job_file
   
   
   subroutine store_star_job_controls(s, ierr)
      use star_private_def
      type (star_info), pointer :: s
      integer, intent(out) :: ierr
      
      ierr = 0
      
      s% job% mesa_dir = mesa_dir
      s% job% eosDT_cache_dir = eosDT_cache_dir
      s% job% kap_cache_dir = kap_cache_dir
      s% job% rates_cache_dir = rates_cache_dir
      s% job% pause_before_terminate = pause_before_terminate
      s% job% profile_columns_file = profile_columns_file
      s% job% history_columns_file = history_columns_file
      s% job% show_log_description_at_start = show_log_description_at_start
      s% job% list_net_reactions = list_net_reactions
      s% job% show_net_reactions_info = show_net_reactions_info
      s% job% show_net_species_info = show_net_species_info
      s% job% first_model_for_timing = first_model_for_timing
      s% job% steps_before_start_timing = steps_before_start_timing
      s% job% show_eqns_and_vars_names = show_eqns_and_vars_names
      s% job% pgstar_flag = pgstar_flag
      s% job% disable_pgstar_during_relax_flag = disable_pgstar_during_relax_flag
      s% job% clear_initial_pgstar_history = clear_initial_pgstar_history
      s% job% clear_pgstar_history = clear_pgstar_history
      s% job% save_pgstar_files_when_terminate = save_pgstar_files_when_terminate
      s% job% save_photo_when_terminate = save_photo_when_terminate
      s% job% load_saved_photo = load_saved_photo
      s% job% saved_photo_name = saved_photo_name
      s% job% save_photo_filename = save_photo_filename
      s% job% save_photo_number = save_photo_number
      s% job% write_profile_when_terminate = write_profile_when_terminate
      s% job% show_retry_counts_when_terminate = show_retry_counts_when_terminate
      s% job% show_timestep_limit_counts_when_terminate = show_timestep_limit_counts_when_terminate
      s% job% filename_for_profile_when_terminate = filename_for_profile_when_terminate
      s% job% create_pre_main_sequence_model = create_pre_main_sequence_model
      s% job% pre_ms_relax_to_start_radiative_core = pre_ms_relax_to_start_radiative_core
      s% job% pre_ms_relax_num_steps = pre_ms_relax_num_steps
      s% job% pre_ms_min_steps_before_check_radiative_core = pre_ms_min_steps_before_check_radiative_core
      s% job% pre_ms_check_radiative_core_start = pre_ms_check_radiative_core_start
      s% job% pre_ms_check_radiative_core_stop = pre_ms_check_radiative_core_stop
      s% job% pre_ms_check_radiative_core_Lnuc_div_L_limit = pre_ms_check_radiative_core_Lnuc_div_L_limit
      s% job% pre_ms_check_radiative_core_min_mass = pre_ms_check_radiative_core_min_mass
      s% job% pre_ms_T_c = pre_ms_T_c
      s% job% pre_ms_guess_rho_c = pre_ms_guess_rho_c
      s% job% pre_ms_d_log10_P = pre_ms_d_log10_P
      s% job% pre_ms_logT_surf_limit = pre_ms_logT_surf_limit
      s% job% pre_ms_logP_surf_limit = pre_ms_logP_surf_limit
      s% job% create_initial_model = create_initial_model
      s% job% initial_model_relax_num_steps = initial_model_relax_num_steps
      s% job% radius_in_cm_for_create_initial_model = radius_in_cm_for_create_initial_model
      s% job% mass_in_gm_for_create_initial_model = mass_in_gm_for_create_initial_model
      
      s% job% center_logP_1st_try_for_create_initial_model = center_logP_1st_try_for_create_initial_model
      s% job% entropy_1st_try_for_create_initial_model = entropy_1st_try_for_create_initial_model
      s% job% abs_e01_tolerance_for_create_initial_model = abs_e01_tolerance_for_create_initial_model
      s% job% abs_e02_tolerance_for_create_initial_model = abs_e02_tolerance_for_create_initial_model
      s% job% max_tries_for_create_initial_model = max_tries_for_create_initial_model
      
      s% job% initial_model_eps = initial_model_eps
      s% job% save_star_job_namelist = save_star_job_namelist
      s% job% star_job_namelist_name = star_job_namelist_name
      s% job% echo_at_start = echo_at_start
      s% job% echo_at_end = echo_at_end
      s% job% load_saved_model = load_saved_model
      s% job% load_model_filename = load_model_filename
      s% job% create_merger_model = create_merger_model
      s% job% saved_model_for_merger_1 = saved_model_for_merger_1
      s% job% saved_model_for_merger_2 = saved_model_for_merger_2
      s% job% set_max_dt_to_frac_lifetime = set_max_dt_to_frac_lifetime
      s% job% max_frac_of_lifetime_per_step = max_frac_of_lifetime_per_step
      s% job% astero_just_call_my_extras_check_model = astero_just_call_my_extras_check_model
      s% job% relax_mass = relax_mass
      s% job% relax_initial_mass = relax_initial_mass
      s% job% new_mass = new_mass
      s% job% lg_max_abs_mdot = lg_max_abs_mdot
      s% job% relax_mass_to_remove_H_env = relax_mass_to_remove_H_env
      s% job% relax_initial_mass_to_remove_H_env = relax_initial_mass_to_remove_H_env
      s% job% extra_mass_retained_by_remove_H_env = extra_mass_retained_by_remove_H_env
      s% job% relax_mass_scale = relax_mass_scale
      s% job% relax_initial_mass_scale = relax_initial_mass_scale
      s% job% dlgm_per_step = dlgm_per_step
      
      s% job% relax_M_center_dt = relax_M_center_dt
      s% job% change_mass_years_for_dt = change_mass_years_for_dt
      s% job% relax_M_center = relax_M_center
      s% job% relax_initial_M_center = relax_initial_M_center
      s% job% relax_core = relax_core
      s% job% relax_initial_core = relax_initial_core
      s% job% new_core_mass = new_core_mass
      s% job% dlg_core_mass_per_step = dlg_core_mass_per_step
      s% job% relax_core_years_for_dt = relax_core_years_for_dt
      s% job% core_avg_rho = core_avg_rho
      s% job% core_avg_eps = core_avg_eps
      
      s% job% relax_R_center = relax_R_center
      s% job% relax_initial_R_center = relax_initial_R_center
      s% job% new_R_center = new_R_center
      s% job% dlgR_per_step = dlgR_per_step
      s% job% relax_R_center_dt = relax_R_center_dt
      
      s% job% set_v_center = set_v_center
      s% job% set_initial_v_center = set_initial_v_center
      
      s% job% relax_v_center = relax_v_center
      s% job% relax_initial_v_center = relax_initial_v_center
      s% job% new_v_center = new_v_center
      s% job% dv_per_step = dv_per_step
      s% job% relax_v_center_dt = relax_v_center_dt
      
      s% job% zero_alpha_RTI = zero_alpha_RTI
      s% job% zero_initial_alpha_RTI = zero_initial_alpha_RTI
      
      s% job% set_L_center = set_L_center
      s% job% set_initial_L_center = set_initial_L_center
      s% job% relax_L_center = relax_L_center
      s% job% relax_initial_L_center = relax_initial_L_center
      s% job% new_L_center = new_L_center
      s% job% dlgL_per_step = dlgL_per_step
      s% job% relax_L_center_dt = relax_L_center_dt
      
      s% job% remove_initial_center_at_cell_k = remove_initial_center_at_cell_k
      s% job% remove_initial_center_by_temperature = remove_initial_center_by_temperature
      s% job% remove_initial_center_by_mass_fraction_q = remove_initial_center_by_mass_fraction_q
      s% job% remove_initial_center_by_mass_gm = remove_initial_center_by_mass_gm
      s% job% remove_initial_center_by_delta_mass_Msun = remove_initial_center_by_delta_mass_Msun
      s% job% remove_initial_center_by_delta_mass_gm = remove_initial_center_by_delta_mass_gm
      s% job% remove_initial_center_by_mass_Msun = remove_initial_center_by_mass_Msun
      s% job% remove_initial_center_by_radius_cm = remove_initial_center_by_radius_cm
      s% job% remove_initial_center_by_radius_Rsun = remove_initial_center_by_radius_Rsun
      s% job% remove_initial_center_by_he4 = remove_initial_center_by_he4
      s% job% remove_initial_center_by_c12_o16 = remove_initial_center_by_c12_o16
      s% job% remove_initial_center_by_si28 = remove_initial_center_by_si28
      s% job% remove_initial_center_to_reduce_co56_ni56 = remove_initial_center_to_reduce_co56_ni56
      s% job% remove_initial_center_by_ye = remove_initial_center_by_ye
      s% job% remove_initial_center_by_entropy = remove_initial_center_by_entropy
      s% job% remove_initial_center_by_infall_kms = remove_initial_center_by_infall_kms
      s% job% remove_initial_center_at_inner_max_abs_v = remove_initial_center_at_inner_max_abs_v
      s% job% remove_initial_fe_core = remove_initial_fe_core
      
      s% job% remove_center_at_cell_k = remove_center_at_cell_k
      s% job% remove_center_by_temperature = remove_center_by_temperature
      s% job% remove_center_by_mass_fraction_q = remove_center_by_mass_fraction_q
      s% job% remove_center_by_delta_mass_gm = remove_center_by_delta_mass_gm
      s% job% remove_center_by_delta_mass_Msun = remove_center_by_delta_mass_Msun
      s% job% remove_center_by_mass_gm = remove_center_by_mass_gm
      s% job% remove_center_by_mass_Msun = remove_center_by_mass_Msun
      s% job% remove_center_by_radius_cm = remove_center_by_radius_cm
      s% job% remove_center_by_radius_Rsun = remove_center_by_radius_Rsun
      s% job% remove_center_by_he4 = remove_center_by_he4
      s% job% remove_center_by_c12_o16 = remove_center_by_c12_o16
      s% job% remove_center_by_si28 = remove_center_by_si28
      s% job% remove_center_to_reduce_co56_ni56 = remove_center_to_reduce_co56_ni56
      s% job% remove_center_by_ye = remove_center_by_ye
      s% job% remove_center_by_entropy = remove_center_by_entropy
      s% job% remove_center_by_infall_kms = remove_center_by_infall_kms
      s% job% remove_center_at_inner_max_abs_v = remove_center_at_inner_max_abs_v
      s% job% remove_fe_core = remove_fe_core
      s% job% remove_fallback_at_each_step = remove_fallback_at_each_step
      s% job% fallback_check_total_energy = fallback_check_total_energy
      s% job% remove_fallback_speed_limit = remove_fallback_speed_limit
      s% job% remove_center_set_zero_v_center = remove_center_set_zero_v_center
      s% job% retain_fallback_at_each_step = retain_fallback_at_each_step
      s% job% limit_center_logP_at_each_step = limit_center_logP_at_each_step
      s% job% remove_center_adjust_L_center = remove_center_adjust_L_center
      s% job% remove_center_logRho_limit = remove_center_logRho_limit
      
      s% job% zero_initial_inner_v_by_mass_Msun = zero_initial_inner_v_by_mass_Msun
      s% job% zero_inner_v_by_mass_Msun = zero_inner_v_by_mass_Msun
      
      s% job% remove_initial_surface_at_cell_k = remove_initial_surface_at_cell_k
      s% job% remove_initial_surface_at_he_core_boundary = remove_initial_surface_at_he_core_boundary
      s% job% remove_initial_surface_by_optical_depth = remove_initial_surface_by_optical_depth
      s% job% remove_initial_surface_by_density = remove_initial_surface_by_density
      s% job% remove_initial_surface_by_pressure = remove_initial_surface_by_pressure
      s% job% remove_initial_surface_by_mass_fraction_q = remove_initial_surface_by_mass_fraction_q
      s% job% remove_initial_surface_by_mass_gm = remove_initial_surface_by_mass_gm
      s% job% remove_initial_surface_by_radius_cm = remove_initial_surface_by_radius_cm
      s% job% remove_initial_surface_by_mass_Msun = remove_initial_surface_by_mass_Msun
      s% job% remove_initial_surface_by_radius_Rsun = remove_initial_surface_by_radius_Rsun
      s% job% remove_initial_surface_by_v_surf_km_s = remove_initial_surface_by_v_surf_km_s
      s% job% remove_initial_surface_by_v_surf_div_cs = remove_initial_surface_by_v_surf_div_cs
      s% job% remove_initial_surface_by_v_surf_div_v_escape = remove_initial_surface_by_v_surf_div_v_escape
      
      s% job% remove_surface_at_cell_k = remove_surface_at_cell_k
      s% job% remove_surface_at_he_core_boundary = remove_surface_at_he_core_boundary
      s% job% remove_surface_by_optical_depth = remove_surface_by_optical_depth
      s% job% remove_surface_by_density = remove_surface_by_density
      s% job% remove_surface_by_pressure = remove_surface_by_pressure
      s% job% remove_surface_by_mass_fraction_q = remove_surface_by_mass_fraction_q
      s% job% remove_surface_by_mass_gm = remove_surface_by_mass_gm
      s% job% remove_surface_by_radius_cm = remove_surface_by_radius_cm
      s% job% remove_surface_by_mass_Msun = remove_surface_by_mass_Msun
      s% job% remove_surface_by_radius_Rsun = remove_surface_by_radius_Rsun
      s% job% remove_surface_by_v_surf_km_s = remove_surface_by_v_surf_km_s
      s% job% remove_surface_by_v_surf_div_cs = remove_surface_by_v_surf_div_cs
      s% job% remove_surface_by_v_surf_div_v_escape = remove_surface_by_v_surf_div_v_escape
      s% job% min_q_for_remove_surface_by_v_surf_div_v_escape = min_q_for_remove_surface_by_v_surf_div_v_escape
      s% job% max_q_for_remove_surface_by_v_surf_div_v_escape = max_q_for_remove_surface_by_v_surf_div_v_escape
      
      s% job% remove_surface_do_jrot = remove_surface_do_jrot
      s% job% remove_surface_do_entropy = remove_surface_do_entropy
      s% job% remove_surface_turn_off_energy_sources_and_sinks = remove_surface_turn_off_energy_sources_and_sinks
      s% job% remove_surface_by_relax_to_star_cut = remove_surface_by_relax_to_star_cut
      
      s% job% report_mass_not_fe56 = report_mass_not_fe56
      s% job% relax_dxdt_nuc_factor = relax_dxdt_nuc_factor
      s% job% relax_initial_dxdt_nuc_factor = relax_initial_dxdt_nuc_factor
      s% job% new_dxdt_nuc_factor = new_dxdt_nuc_factor
      s% job% dxdt_nuc_factor_multiplier = dxdt_nuc_factor_multiplier
      s% job% relax_eps_nuc_factor = relax_eps_nuc_factor
      s% job% relax_initial_eps_nuc_factor = relax_initial_eps_nuc_factor
      s% job% new_eps_nuc_factor = new_eps_nuc_factor
      s% job% eps_nuc_factor_multiplier = eps_nuc_factor_multiplier
      s% job% relax_opacity_max = relax_opacity_max
      s% job% relax_initial_opacity_max = relax_initial_opacity_max
      s% job% new_opacity_max = new_opacity_max
      s% job% opacity_max_multiplier = opacity_max_multiplier
      s% job% relax_max_surf_dq = relax_max_surf_dq
      s% job% relax_initial_max_surf_dq = relax_initial_max_surf_dq
      s% job% new_max_surf_dq = new_max_surf_dq
      s% job% max_surf_dq_multiplier = max_surf_dq_multiplier
      
      s% job% relax_tau_factor = relax_tau_factor
      s% job% relax_initial_tau_factor = relax_initial_tau_factor
      s% job% set_tau_factor = set_tau_factor
      s% job% set_initial_tau_factor = set_initial_tau_factor
      s% job% relax_to_this_tau_factor = relax_to_this_tau_factor
      s% job% set_to_this_tau_factor = set_to_this_tau_factor
      s% job% dlogtau_factor = dlogtau_factor
      s% job% set_tau_factor_after_core_He_burn = set_tau_factor_after_core_He_burn
      s% job% set_tau_factor_after_core_C_burn = set_tau_factor_after_core_C_burn
      s% job% relax_tau_factor_after_core_He_burn = relax_tau_factor_after_core_He_burn
      s% job% relax_tau_factor_after_core_C_burn = relax_tau_factor_after_core_C_burn
      
      s% job% adjust_tau_factor_to_surf_density = adjust_tau_factor_to_surf_density
      s% job% base_for_adjust_tau_factor_to_surf_density = base_for_adjust_tau_factor_to_surf_density
      
      s% job% relax_to_this_opacity_factor = relax_to_this_opacity_factor
      s% job% d_opacity_factor = d_opacity_factor
      s% job% relax_opacity_factor = relax_opacity_factor
      s% job% relax_initial_opacity_factor = relax_initial_opacity_factor
      
      s% job% relax_Tsurf_factor = relax_Tsurf_factor
      s% job% relax_initial_Tsurf_factor = relax_initial_Tsurf_factor
      s% job% set_Tsurf_factor = set_Tsurf_factor
      s% job% set_initial_Tsurf_factor = set_initial_Tsurf_factor
      s% job% relax_to_this_Tsurf_factor = relax_to_this_Tsurf_factor
      s% job% set_to_this_Tsurf_factor = set_to_this_Tsurf_factor
      s% job% dlogTsurf_factor = dlogTsurf_factor
      
      s% job% relax_irradiation = relax_irradiation
      s% job% relax_initial_irradiation = relax_initial_irradiation
      s% job% set_irradiation = set_irradiation
      s% job% set_initial_irradiation = set_initial_irradiation
      s% job% relax_irradiation_min_steps = relax_irradiation_min_steps
      s% job% relax_to_this_irrad_flux = relax_to_this_irrad_flux
      s% job% set_to_this_irrad_flux = set_to_this_irrad_flux
      s% job% irrad_col_depth = irrad_col_depth
      s% job% relax_irradiation_max_yrs_dt = relax_irradiation_max_yrs_dt
      s% job% relax_mass_change = relax_mass_change
      s% job% relax_initial_mass_change = relax_initial_mass_change
      s% job% relax_mass_change_min_steps = relax_mass_change_min_steps
      s% job% relax_mass_change_max_yrs_dt = relax_mass_change_max_yrs_dt
      s% job% relax_mass_change_init_mdot = relax_mass_change_init_mdot
      s% job% relax_mass_change_final_mdot = relax_mass_change_final_mdot
      s% job% change_RTI_flag = change_RTI_flag
      s% job% change_initial_RTI_flag = change_initial_RTI_flag
      s% job% new_RTI_flag = new_RTI_flag
      s% job% change_RSP_flag = change_RSP_flag
      s% job% change_initial_RSP_flag = change_initial_RSP_flag
      s% job% new_RSP_flag = new_RSP_flag
      s% job% change_RSP2_flag = change_RSP2_flag
      s% job% change_initial_RSP2_flag = change_initial_RSP2_flag
      s% job% change_RSP2_flag_at_model_number = change_RSP2_flag_at_model_number
      s% job% new_RSP2_flag = new_RSP2_flag
      s% job% create_RSP2_model = create_RSP2_model
      s% job% change_w_div_wc_flag = change_w_div_wc_flag
      s% job% change_initial_w_div_wc_flag = change_initial_w_div_wc_flag
      s% job% new_w_div_wc_flag = new_w_div_wc_flag
      s% job% change_j_rot_flag = change_j_rot_flag
      s% job% change_initial_j_rot_flag = change_initial_j_rot_flag
      s% job% new_j_rot_flag = new_j_rot_flag
      
      s% job% create_RSP_model = create_RSP_model
      
      s% job% change_v_flag = change_v_flag
      s% job% change_initial_v_flag = change_initial_v_flag
      s% job% new_v_flag = new_v_flag
      s% job% change_D_omega_flag = change_D_omega_flag
      s% job% change_initial_D_omega_flag = change_initial_D_omega_flag
      s% job% new_D_omega_flag = new_D_omega_flag
      s% job% change_am_nu_rot_flag = change_am_nu_rot_flag
      s% job% change_initial_am_nu_rot_flag = change_initial_am_nu_rot_flag
      s% job% new_am_nu_rot_flag = new_am_nu_rot_flag
      s% job% use_D_omega_for_am_nu_rot = use_D_omega_for_am_nu_rot
      
      s% job% change_u_flag = change_u_flag
      s% job% change_initial_u_flag = change_initial_u_flag
      s% job% new_u_flag = new_u_flag
      
      s% job% change_reconstruction_flag = change_reconstruction_flag
      s% job% change_initial_reconstruction_flag = change_initial_reconstruction_flag
      s% job% new_reconstruction_flag = new_reconstruction_flag
      
      s% job% center_ye_limit_for_v_flag = center_ye_limit_for_v_flag
      s% job% change_rotation_flag = change_rotation_flag
      s% job% change_initial_rotation_flag = change_initial_rotation_flag
      s% job% new_rotation_flag = new_rotation_flag
      s% job% use_w_div_wc_flag_with_rotation = use_w_div_wc_flag_with_rotation
      s% job% use_j_rot_flag_with_rotation = use_j_rot_flag_with_rotation
      s% job% set_omega = set_omega
      s% job% set_initial_omega = set_initial_omega
      s% job% set_omega_step_limit = set_omega_step_limit
      s% job% set_near_zams_omega_steps = set_near_zams_omega_steps
      s% job% new_omega = new_omega
      s% job% set_omega_div_omega_crit = set_omega_div_omega_crit
      s% job% set_initial_omega_div_omega_crit = set_initial_omega_div_omega_crit
      s% job% set_omega_div_omega_crit_step_limit = set_omega_div_omega_crit_step_limit
      s% job% set_near_zams_omega_div_omega_crit_steps = set_near_zams_omega_div_omega_crit_steps
      s% job% new_omega_div_omega_crit = new_omega_div_omega_crit
      s% job% set_surface_rotation_v = set_surface_rotation_v
      s% job% set_initial_surface_rotation_v = set_initial_surface_rotation_v
      s% job% set_surf_rotation_v_step_limit = set_surf_rotation_v_step_limit
      s% job% set_near_zams_surface_rotation_v_steps = set_near_zams_surface_rotation_v_steps
      s% job% new_surface_rotation_v = new_surface_rotation_v
      s% job% relax_omega = relax_omega
      s% job% relax_initial_omega = relax_initial_omega
      s% job% near_zams_relax_omega = near_zams_relax_omega
      s% job% relax_omega_div_omega_crit = relax_omega_div_omega_crit
      s% job% relax_initial_omega_div_omega_crit = relax_initial_omega_div_omega_crit
      s% job% near_zams_relax_omega_div_omega_crit = near_zams_relax_omega_div_omega_crit
      s% job% relax_surface_rotation_v = relax_surface_rotation_v
      s% job% relax_initial_surface_rotation_v = relax_initial_surface_rotation_v
      s% job% near_zams_relax_initial_surface_rotation_v = near_zams_relax_initial_surface_rotation_v
      s% job% num_steps_to_relax_rotation = num_steps_to_relax_rotation
      s% job% relax_omega_max_yrs_dt = relax_omega_max_yrs_dt
      s% job% set_uniform_initial_composition = set_uniform_initial_composition
      s% job% initial_h1 = initial_h1
      s% job% initial_h2 = initial_h2
      s% job% initial_he3 = initial_he3
      s% job% initial_he4 = initial_he4
      s% job% initial_zfracs = initial_zfracs
      s% job% dump_missing_metals_into_heaviest = dump_missing_metals_into_heaviest
      s% job% relax_initial_composition = relax_initial_composition
      s% job% relax_initial_to_xaccrete = relax_initial_to_xaccrete
      s% job% relax_composition_filename = relax_composition_filename
      s% job% num_steps_to_relax_composition = num_steps_to_relax_composition
      s% job% timescale_for_relax_composition = timescale_for_relax_composition
      s% job% relax_initial_angular_momentum = relax_initial_angular_momentum
      s% job% max_steps_to_relax_angular_momentum = max_steps_to_relax_angular_momentum
      s% job% timescale_for_relax_angular_momentum = timescale_for_relax_angular_momentum
      s% job% max_dt_for_relax_angular_momentum = max_dt_for_relax_angular_momentum
      s% job% num_timescales_for_relax_angular_momentum = num_timescales_for_relax_angular_momentum
      s% job% relax_angular_momentum_filename = relax_angular_momentum_filename
      s% job% relax_angular_momentum_constant_omega_center = relax_angular_momentum_constant_omega_center
      s% job% relax_initial_entropy = relax_initial_entropy
      s% job% max_steps_to_relax_entropy = max_steps_to_relax_entropy
      s% job% timescale_for_relax_entropy = timescale_for_relax_entropy
      s% job% max_dt_for_relax_entropy = max_dt_for_relax_entropy
      s% job% num_timescales_for_relax_entropy = num_timescales_for_relax_entropy
      s% job% relax_entropy_filename = relax_entropy_filename
      s% job% get_entropy_for_relax_from_eos = get_entropy_for_relax_from_eos
      s% job% report_cell_for_xm = report_cell_for_xm
      s% job% set_to_xa_for_accretion = set_to_xa_for_accretion
      s% job% set_initial_to_xa_for_accretion = set_initial_to_xa_for_accretion
      s% job% set_nzlo = set_nzlo
      s% job% set_nzhi = set_nzhi
      s% job% change_Y = change_Y
      s% job% change_initial_Y = change_initial_Y
      s% job% relax_Y = relax_Y
      s% job% relax_initial_Y = relax_initial_Y
      s% job% relax_Y_minq = relax_Y_minq
      s% job% relax_Y_maxq = relax_Y_maxq
      s% job% new_Y = new_Y
      s% job% change_Z = change_Z
      s% job% change_initial_Z = change_initial_Z
      s% job% relax_Z = relax_Z
      s% job% relax_initial_Z = relax_initial_Z
      s% job% relax_Z_minq = relax_Z_minq
      s% job% relax_Z_maxq = relax_Z_maxq
      s% job% new_Z = new_Z
      s% job% steps_to_take_before_terminate = steps_to_take_before_terminate
      s% job% stop_if_this_file_exists = stop_if_this_file_exists
      s% job% set_initial_age = set_initial_age
      s% job% initial_age = initial_age
      s% job% set_initial_model_number = set_initial_model_number
      s% job% initial_model_number = initial_model_number
      s% job% set_initial_number_retries = set_initial_number_retries
      s% job% initial_number_retries = initial_number_retries
      s% job% set_initial_dt = set_initial_dt
      s% job% limit_initial_dt = limit_initial_dt
      s% job% years_for_initial_dt = years_for_initial_dt
      s% job% seconds_for_initial_dt = seconds_for_initial_dt
      
      s% job% set_initial_cumulative_energy_error = set_initial_cumulative_energy_error
      s% job% set_cumulative_energy_error = set_cumulative_energy_error
      s% job% set_cumulative_energy_error_at_step = set_cumulative_energy_error_at_step
      s% job% set_cumulative_energy_error_each_step_if_age_less_than = set_cumulative_energy_error_each_step_if_age_less_than
      s% job% new_cumulative_energy_error = new_cumulative_energy_error
      s% job% set_cumulative_energy_error_each_relax = set_cumulative_energy_error_each_relax
      
      s% job% change_net = change_net
      s% job% change_initial_net = change_initial_net
      s% job% new_net_name = new_net_name
      s% job% change_small_net = change_small_net
      s% job% change_initial_small_net = change_initial_small_net
      s% job% new_small_net_name = new_small_net_name
      
      s% job% h_he_net = h_he_net
      s% job% co_net = co_net
      s% job% adv_net = adv_net
      s% job% adjust_abundances_for_new_isos = adjust_abundances_for_new_isos
      s% job% set_uniform_xa_from_file = set_uniform_xa_from_file
      s% job% set_uniform_initial_xa_from_file = set_uniform_initial_xa_from_file
      s% job% file_for_uniform_xa = file_for_uniform_xa
      
      s% job% mix_section = mix_section
      s% job% mix_initial_section = mix_initial_section
      s% job% mix_section_nzlo = mix_section_nzlo
      s% job% mix_section_nzhi = mix_section_nzhi
      
      s% job% T9_weaklib_full_off = T9_weaklib_full_off
      s% job% T9_weaklib_full_on = T9_weaklib_full_on
      s% job% weaklib_blend_hi_Z = weaklib_blend_hi_Z
      s% job% T9_weaklib_full_off_hi_Z = T9_weaklib_full_off_hi_Z
      s% job% T9_weaklib_full_on_hi_Z = T9_weaklib_full_on_hi_Z
      
      s% job% use_suzuki_weak_rates = use_suzuki_weak_rates
      s% job% use_3a_fl87 = use_3a_fl87
      
      s% job% use_special_weak_rates = use_special_weak_rates
      s% job% special_weak_states_file = special_weak_states_file
      s% job% special_weak_transitions_file = special_weak_transitions_file
      s% job% ion_coulomb_corrections = ion_coulomb_corrections
      s% job% electron_coulomb_corrections = electron_coulomb_corrections
      
      s% job% mix_envelope_down_to_T = mix_envelope_down_to_T
      s% job% mix_initial_envelope_down_to_T = mix_initial_envelope_down_to_T
      s% job% auto_extend_net = auto_extend_net
      
      s% job% enable_adaptive_network = enable_adaptive_network
      s% job% min_x_for_keep = min_x_for_keep
      s% job% min_x_for_n = min_x_for_n
      s% job% min_x_for_add = min_x_for_add
      s% job% max_Z_for_add = max_Z_for_add
      s% job% max_N_for_add = max_N_for_add
      s% job% max_A_for_add = max_A_for_add
      
      s% job% save_model_number = save_model_number
      s% job% save_model_filename = save_model_filename
      s% job% save_model_when_terminate = save_model_when_terminate
      s% job% required_termination_code_string = required_termination_code_string
      s% job% profile_starting_model = profile_starting_model
      s% job% profile_model_number = profile_model_number
      s% job% report_retries = report_retries
      
      s% job% net_reaction_filename = net_reaction_filename
      s% job% jina_reaclib_filename = jina_reaclib_filename
      s% job% jina_reaclib_min_T9 = jina_reaclib_min_T9
      s% job% rate_tables_dir = rate_tables_dir
      s% job% rate_cache_suffix = rate_cache_suffix
      s% job% read_extra_star_job_inlist1 = read_extra_star_job_inlist1
      s% job% extra_star_job_inlist1_name = extra_star_job_inlist1_name
      s% job% read_extra_star_job_inlist2 = read_extra_star_job_inlist2
      s% job% extra_star_job_inlist2_name = extra_star_job_inlist2_name
      s% job% read_extra_star_job_inlist3 = read_extra_star_job_inlist3
      s% job% extra_star_job_inlist3_name = extra_star_job_inlist3_name
      s% job% read_extra_star_job_inlist4 = read_extra_star_job_inlist4
      s% job% extra_star_job_inlist4_name = extra_star_job_inlist4_name
      s% job% read_extra_star_job_inlist5 = read_extra_star_job_inlist5
      s% job% extra_star_job_inlist5_name = extra_star_job_inlist5_name
      s% job% set_abundance_nzlo = set_abundance_nzlo
      s% job% set_abundance_nzhi = set_abundance_nzhi
      s% job% set_abundance = set_abundance
      s% job% set_initial_abundance = set_initial_abundance
      s% job% chem_name = chem_name
      s% job% new_frac = new_frac
      s% job% set_abundance_nzlo = set_abundance_nzlo
      s% job% set_abundance_nzhi = set_abundance_nzhi
      s% job% replace_element = replace_element
      s% job% replace_initial_element = replace_initial_element
      s% job% chem_name1 = chem_name1
      s% job% chem_name2 = chem_name2
      s% job% replace_element_nzlo = replace_element_nzlo
      s% job% replace_element_nzhi = replace_element_nzhi
      s% job% do_special_test = do_special_test
      
      s% job% save_pulse_data_for_model_number = save_pulse_data_for_model_number
      s% job% save_pulse_data_when_terminate = save_pulse_data_when_terminate
      s% job% save_pulse_data_filename = save_pulse_data_filename
      
      s% job% chem_isotopes_filename = chem_isotopes_filename
      s% job% extras_lipar = extras_lipar
      s% job% extras_lrpar = extras_lrpar
      s% job% extras_lcpar = extras_lcpar
      s% job% extras_llpar = extras_llpar
      s% job% extras_ipar = extras_ipar
      s% job% extras_rpar = extras_rpar
      s% job% extras_cpar = extras_cpar
      s% job% extras_lpar = extras_lpar
      s% job% num_special_rate_factors = num_special_rate_factors
      s% job% special_rate_factor = special_rate_factor
      s% job% filename_of_special_rate = filename_of_special_rate
      s% job% reaction_for_special_factor = reaction_for_special_factor
      s% job% color_num_files = color_num_files
      s% job% color_file_names = color_file_names
      s% job% color_num_colors = color_num_colors
      
      s% job% warn_run_star_extras = warn_run_star_extras
      s% job% report_garbage_collection = report_garbage_collection
      s% job% num_steps_for_garbage_collection = num_steps_for_garbage_collection
   
   end subroutine store_star_job_controls
   
   
   subroutine set_default_star_job_controls
      required_termination_code_string(:) = ''
      extras_ipar(:) = 0
      extras_rpar(:) = 0
      extras_cpar(:) = ''
      extras_lpar(:) = .false.
      special_rate_factor(:) = 1d0
      filename_of_special_rate(:) = ''
      reaction_for_special_factor(:) = ''
      color_num_colors(:) = 0
      color_file_names(:) = ''
      include 'star_job.defaults'
      include 'star_job_dev.defaults'
   end subroutine set_default_star_job_controls
   
   
   subroutine check_star_job_controls(s, ierr)
      type (star_info), pointer :: s
      integer, intent(out) :: ierr
      
      ierr = 0
      
      if (s% job% write_profile_when_terminate) then
         if (len_trim(s% job% filename_for_profile_when_terminate) == 0) then
            write(*, *) "when write_profile_when_terminate = .true.,"
            write(*, *) "filename_for_profile_when_terminate must be non empty"
            ierr = -1
            return
         end if
      end if
      
      if (s% job% save_model_when_terminate) then
         if (len_trim(s% job% save_model_filename) == 0) then
            write(*, *) "when save_model_when_terminate = .true.,"
            write(*, *) "filename_for_profile_when_terminate must be non empty"
            ierr = -1
            return
         end if
      end if
   
   end subroutine check_star_job_controls
   
   
   subroutine set_star_job_controls_for_writing(s, ierr)
      type (star_info), pointer :: s
      integer, intent(out) :: ierr
      
      ierr = 0
      
      mesa_dir = s% job% mesa_dir
      eosDT_cache_dir = s% job% eosDT_cache_dir
      kap_cache_dir = s% job% kap_cache_dir
      rates_cache_dir = s% job% rates_cache_dir
      pause_before_terminate = s% job% pause_before_terminate
      profile_columns_file = s% job% profile_columns_file
      history_columns_file = s% job% history_columns_file
      show_log_description_at_start = s% job% show_log_description_at_start
      list_net_reactions = s% job% list_net_reactions
      show_net_reactions_info = s% job% show_net_reactions_info
      show_net_species_info = s% job% show_net_species_info
      first_model_for_timing = s% job% first_model_for_timing
      steps_before_start_timing = s% job% steps_before_start_timing
      show_eqns_and_vars_names = s% job% show_eqns_and_vars_names
      pgstar_flag = s% job% pgstar_flag
      disable_pgstar_during_relax_flag = s% job% disable_pgstar_during_relax_flag
      clear_initial_pgstar_history = s% job% clear_initial_pgstar_history
      clear_pgstar_history = s% job% clear_pgstar_history
      save_pgstar_files_when_terminate = s% job% save_pgstar_files_when_terminate
      save_photo_when_terminate = s% job% save_photo_when_terminate
      load_saved_photo = s% job% load_saved_photo
      saved_photo_name = s% job% saved_photo_name
      save_photo_filename = s% job% save_photo_filename
      save_photo_number = s% job% save_photo_number
      write_profile_when_terminate = s% job% write_profile_when_terminate
      show_retry_counts_when_terminate = s% job% show_retry_counts_when_terminate
      show_timestep_limit_counts_when_terminate = s% job% show_timestep_limit_counts_when_terminate
      filename_for_profile_when_terminate = s% job% filename_for_profile_when_terminate
      create_pre_main_sequence_model = s% job% create_pre_main_sequence_model
      pre_ms_relax_to_start_radiative_core = s% job% pre_ms_relax_to_start_radiative_core
      pre_ms_relax_num_steps = s% job% pre_ms_relax_num_steps
      pre_ms_min_steps_before_check_radiative_core = s% job% pre_ms_min_steps_before_check_radiative_core
      pre_ms_check_radiative_core_start = s% job% pre_ms_check_radiative_core_start
      pre_ms_check_radiative_core_stop = s% job% pre_ms_check_radiative_core_stop
      pre_ms_check_radiative_core_Lnuc_div_L_limit = s% job% pre_ms_check_radiative_core_Lnuc_div_L_limit
      pre_ms_check_radiative_core_min_mass = s% job% pre_ms_check_radiative_core_min_mass
      pre_ms_T_c = s% job% pre_ms_T_c
      pre_ms_guess_rho_c = s% job% pre_ms_guess_rho_c
      pre_ms_d_log10_P = s% job% pre_ms_d_log10_P
      pre_ms_logT_surf_limit = s% job% pre_ms_logT_surf_limit
      pre_ms_logP_surf_limit = s% job% pre_ms_logP_surf_limit
      create_initial_model = s% job% create_initial_model
      initial_model_relax_num_steps = s% job% initial_model_relax_num_steps
      radius_in_cm_for_create_initial_model = s% job% radius_in_cm_for_create_initial_model
      mass_in_gm_for_create_initial_model = s% job% mass_in_gm_for_create_initial_model
      
      center_logP_1st_try_for_create_initial_model = s% job% center_logP_1st_try_for_create_initial_model
      entropy_1st_try_for_create_initial_model = s% job% entropy_1st_try_for_create_initial_model
      abs_e01_tolerance_for_create_initial_model = s% job% abs_e01_tolerance_for_create_initial_model
      abs_e02_tolerance_for_create_initial_model = s% job% abs_e02_tolerance_for_create_initial_model
      max_tries_for_create_initial_model = s% job% max_tries_for_create_initial_model
      
      initial_model_eps = s% job% initial_model_eps
      save_star_job_namelist = s% job% save_star_job_namelist
      star_job_namelist_name = s% job% star_job_namelist_name
      echo_at_start = s% job% echo_at_start
      echo_at_end = s% job% echo_at_end
      load_saved_model = s% job% load_saved_model
      load_model_filename = s% job% load_model_filename
      create_merger_model = s% job% create_merger_model
      saved_model_for_merger_1 = s% job% saved_model_for_merger_1
      saved_model_for_merger_2 = s% job% saved_model_for_merger_2
      set_max_dt_to_frac_lifetime = s% job% set_max_dt_to_frac_lifetime
      max_frac_of_lifetime_per_step = s% job% max_frac_of_lifetime_per_step
      astero_just_call_my_extras_check_model = s% job% astero_just_call_my_extras_check_model
      relax_mass = s% job% relax_mass
      relax_initial_mass = s% job% relax_initial_mass
      new_mass = s% job% new_mass
      lg_max_abs_mdot = s% job% lg_max_abs_mdot
      relax_mass_to_remove_H_env = s% job% relax_mass_to_remove_H_env
      relax_initial_mass_to_remove_H_env = s% job% relax_initial_mass_to_remove_H_env
      extra_mass_retained_by_remove_H_env = s% job% extra_mass_retained_by_remove_H_env
      relax_mass_scale = s% job% relax_mass_scale
      relax_initial_mass_scale = s% job% relax_initial_mass_scale
      dlgm_per_step = s% job% dlgm_per_step
      
      relax_M_center_dt = s% job% relax_M_center_dt
      change_mass_years_for_dt = s% job% change_mass_years_for_dt
      relax_M_center = s% job% relax_M_center
      relax_initial_M_center = s% job% relax_initial_M_center
      relax_core = s% job% relax_core
      relax_initial_core = s% job% relax_initial_core
      new_core_mass = s% job% new_core_mass
      dlg_core_mass_per_step = s% job% dlg_core_mass_per_step
      relax_core_years_for_dt = s% job% relax_core_years_for_dt
      core_avg_rho = s% job% core_avg_rho
      core_avg_eps = s% job% core_avg_eps
      
      relax_R_center = s% job% relax_R_center
      relax_initial_R_center = s% job% relax_initial_R_center
      new_R_center = s% job% new_R_center
      dlgR_per_step = s% job% dlgR_per_step
      relax_R_center_dt = s% job% relax_R_center_dt
      
      set_v_center = s% job% set_v_center
      set_initial_v_center = s% job% set_initial_v_center
      
      relax_v_center = s% job% relax_v_center
      relax_initial_v_center = s% job% relax_initial_v_center
      new_v_center = s% job% new_v_center
      dv_per_step = s% job% dv_per_step
      relax_v_center_dt = s% job% relax_v_center_dt
      
      zero_alpha_RTI = s% job% zero_alpha_RTI
      zero_initial_alpha_RTI = s% job% zero_initial_alpha_RTI
      
      set_L_center = s% job% set_L_center
      set_initial_L_center = s% job% set_initial_L_center
      relax_L_center = s% job% relax_L_center
      relax_initial_L_center = s% job% relax_initial_L_center
      new_L_center = s% job% new_L_center
      dlgL_per_step = s% job% dlgL_per_step
      relax_L_center_dt = s% job% relax_L_center_dt
      
      remove_center_at_cell_k = s% job% remove_center_at_cell_k
      remove_center_by_temperature = s% job% remove_center_by_temperature
      remove_center_by_mass_fraction_q = s% job% remove_center_by_mass_fraction_q
      remove_center_by_delta_mass_gm = s% job% remove_center_by_delta_mass_gm
      remove_center_by_delta_mass_Msun = s% job% remove_center_by_delta_mass_Msun
      remove_center_by_mass_gm = s% job% remove_center_by_mass_gm
      remove_center_by_mass_Msun = s% job% remove_center_by_mass_Msun
      remove_center_by_radius_Rsun = s% job% remove_center_by_radius_Rsun
      remove_center_by_radius_cm = s% job% remove_center_by_radius_cm
      remove_center_by_he4 = s% job% remove_center_by_he4
      remove_center_by_c12_o16 = s% job% remove_center_by_c12_o16
      remove_center_by_si28 = s% job% remove_center_by_si28
      remove_center_to_reduce_co56_ni56 = s% job% remove_center_to_reduce_co56_ni56
      remove_center_by_ye = s% job% remove_center_by_ye
      remove_center_by_entropy = s% job% remove_center_by_entropy
      remove_center_by_infall_kms = s% job% remove_center_by_infall_kms
      remove_center_at_inner_max_abs_v = s% job% remove_center_at_inner_max_abs_v
      remove_fe_core = s% job% remove_fe_core
      remove_fallback_at_each_step = s% job% remove_fallback_at_each_step
      fallback_check_total_energy = s% job% fallback_check_total_energy
      remove_fallback_speed_limit = s% job% remove_fallback_speed_limit
      remove_center_set_zero_v_center = s% job% remove_center_set_zero_v_center
      limit_center_logP_at_each_step = s% job% limit_center_logP_at_each_step
      remove_center_adjust_L_center = s% job% remove_center_adjust_L_center
      remove_center_logRho_limit = s% job% remove_center_logRho_limit
      
      remove_initial_center_at_cell_k = s% job% remove_initial_center_at_cell_k
      remove_initial_center_by_temperature = s% job% remove_initial_center_by_temperature
      remove_initial_center_by_mass_fraction_q = s% job% remove_initial_center_by_mass_fraction_q
      remove_initial_center_by_delta_mass_gm = s% job% remove_initial_center_by_delta_mass_gm
      remove_initial_center_by_delta_mass_Msun = s% job% remove_initial_center_by_delta_mass_Msun
      remove_initial_center_by_mass_gm = s% job% remove_initial_center_by_mass_gm
      remove_initial_center_by_mass_Msun = s% job% remove_initial_center_by_mass_Msun
      remove_initial_center_by_radius_Rsun = s% job% remove_initial_center_by_radius_Rsun
      remove_initial_center_by_radius_cm = s% job% remove_initial_center_by_radius_cm
      remove_initial_center_by_he4 = s% job% remove_initial_center_by_he4
      remove_initial_center_by_c12_o16 = s% job% remove_initial_center_by_c12_o16
      remove_initial_center_by_si28 = s% job% remove_initial_center_by_si28
      remove_initial_center_to_reduce_co56_ni56 = s% job% remove_initial_center_to_reduce_co56_ni56
      remove_initial_center_by_ye = s% job% remove_initial_center_by_ye
      remove_initial_center_by_entropy = s% job% remove_initial_center_by_entropy
      remove_initial_center_by_infall_kms = s% job% remove_initial_center_by_infall_kms
      remove_initial_center_at_inner_max_abs_v = s% job% remove_initial_center_at_inner_max_abs_v
      remove_initial_fe_core = s% job% remove_initial_fe_core
      
      zero_initial_inner_v_by_mass_Msun = s% job% zero_initial_inner_v_by_mass_Msun
      zero_inner_v_by_mass_Msun = s% job% zero_inner_v_by_mass_Msun
      
      remove_surface_at_cell_k = s% job% remove_surface_at_cell_k
      remove_surface_at_he_core_boundary = s% job% remove_surface_at_he_core_boundary
      remove_surface_by_optical_depth = s% job% remove_surface_by_optical_depth
      remove_surface_by_density = s% job% remove_surface_by_density
      remove_surface_by_pressure = s% job% remove_surface_by_pressure
      remove_surface_by_mass_fraction_q = s% job% remove_surface_by_mass_fraction_q
      remove_surface_by_mass_gm = s% job% remove_surface_by_mass_gm
      remove_surface_by_radius_cm = s% job% remove_surface_by_radius_cm
      remove_surface_by_mass_Msun = s% job% remove_surface_by_mass_Msun
      remove_surface_by_radius_Rsun = s% job% remove_surface_by_radius_Rsun
      remove_surface_by_v_surf_km_s = s% job% remove_surface_by_v_surf_km_s
      remove_surface_by_v_surf_div_cs = s% job% remove_surface_by_v_surf_div_cs
      remove_surface_by_v_surf_div_v_escape = s% job% remove_surface_by_v_surf_div_v_escape
      min_q_for_remove_surface_by_v_surf_div_v_escape = s% job% min_q_for_remove_surface_by_v_surf_div_v_escape
      max_q_for_remove_surface_by_v_surf_div_v_escape = s% job% max_q_for_remove_surface_by_v_surf_div_v_escape
      
      remove_surface_do_jrot = s% job% remove_surface_do_jrot
      remove_surface_do_entropy = s% job% remove_surface_do_entropy
      remove_surface_turn_off_energy_sources_and_sinks = s% job% remove_surface_turn_off_energy_sources_and_sinks
      remove_surface_by_relax_to_star_cut = s% job% remove_surface_by_relax_to_star_cut
      
      remove_initial_surface_at_cell_k = s% job% remove_initial_surface_at_cell_k
      remove_initial_surface_at_he_core_boundary = s% job% remove_initial_surface_at_he_core_boundary
      remove_initial_surface_by_optical_depth = s% job% remove_initial_surface_by_optical_depth
      remove_initial_surface_by_density = s% job% remove_initial_surface_by_density
      remove_initial_surface_by_pressure = s% job% remove_initial_surface_by_pressure
      remove_initial_surface_by_mass_fraction_q = s% job% remove_initial_surface_by_mass_fraction_q
      remove_initial_surface_by_mass_gm = s% job% remove_initial_surface_by_mass_gm
      remove_initial_surface_by_radius_cm = s% job% remove_initial_surface_by_radius_cm
      remove_initial_surface_by_mass_Msun = s% job% remove_initial_surface_by_mass_Msun
      remove_initial_surface_by_radius_Rsun = s% job% remove_initial_surface_by_radius_Rsun
      remove_initial_surface_by_v_surf_km_s = s% job% remove_initial_surface_by_v_surf_km_s
      remove_initial_surface_by_v_surf_div_cs = s% job% remove_initial_surface_by_v_surf_div_cs
      remove_initial_surface_by_v_surf_div_v_escape = s% job% remove_initial_surface_by_v_surf_div_v_escape
      
      report_mass_not_fe56 = s% job% report_mass_not_fe56
      relax_dxdt_nuc_factor = s% job% relax_dxdt_nuc_factor
      relax_initial_dxdt_nuc_factor = s% job% relax_initial_dxdt_nuc_factor
      new_dxdt_nuc_factor = s% job% new_dxdt_nuc_factor
      dxdt_nuc_factor_multiplier = s% job% dxdt_nuc_factor_multiplier
      relax_eps_nuc_factor = s% job% relax_eps_nuc_factor
      relax_initial_eps_nuc_factor = s% job% relax_initial_eps_nuc_factor
      new_eps_nuc_factor = s% job% new_eps_nuc_factor
      eps_nuc_factor_multiplier = s% job% eps_nuc_factor_multiplier
      relax_opacity_max = s% job% relax_opacity_max
      relax_initial_opacity_max = s% job% relax_initial_opacity_max
      new_opacity_max = s% job% new_opacity_max
      opacity_max_multiplier = s% job% opacity_max_multiplier
      relax_max_surf_dq = s% job% relax_max_surf_dq
      relax_initial_max_surf_dq = s% job% relax_initial_max_surf_dq
      new_max_surf_dq = s% job% new_max_surf_dq
      max_surf_dq_multiplier = s% job% max_surf_dq_multiplier
      
      relax_tau_factor = s% job% relax_tau_factor
      relax_initial_tau_factor = s% job% relax_initial_tau_factor
      set_tau_factor = s% job% set_tau_factor
      set_initial_tau_factor = s% job% set_initial_tau_factor
      relax_to_this_tau_factor = s% job% relax_to_this_tau_factor
      set_to_this_tau_factor = s% job% set_to_this_tau_factor
      dlogtau_factor = s% job% dlogtau_factor
      set_tau_factor_after_core_He_burn = s% job% set_tau_factor_after_core_He_burn
      set_tau_factor_after_core_C_burn = s% job% set_tau_factor_after_core_C_burn
      relax_tau_factor_after_core_He_burn = s% job% relax_tau_factor_after_core_He_burn
      relax_tau_factor_after_core_C_burn = s% job% relax_tau_factor_after_core_C_burn
      
      adjust_tau_factor_to_surf_density = s% job% adjust_tau_factor_to_surf_density
      base_for_adjust_tau_factor_to_surf_density = s% job% base_for_adjust_tau_factor_to_surf_density
      
      relax_to_this_opacity_factor = s% job% relax_to_this_opacity_factor
      d_opacity_factor = s% job% d_opacity_factor
      relax_opacity_factor = s% job% relax_opacity_factor
      relax_initial_opacity_factor = s% job% relax_initial_opacity_factor
      
      relax_Tsurf_factor = s% job% relax_Tsurf_factor
      relax_initial_Tsurf_factor = s% job% relax_initial_Tsurf_factor
      set_Tsurf_factor = s% job% set_Tsurf_factor
      set_initial_Tsurf_factor = s% job% set_initial_Tsurf_factor
      relax_to_this_Tsurf_factor = s% job% relax_to_this_Tsurf_factor
      set_to_this_Tsurf_factor = s% job% set_to_this_Tsurf_factor
      dlogTsurf_factor = s% job% dlogTsurf_factor
      
      relax_irradiation = s% job% relax_irradiation
      relax_initial_irradiation = s% job% relax_initial_irradiation
      set_irradiation = s% job% set_irradiation
      set_initial_irradiation = s% job% set_initial_irradiation
      relax_irradiation_min_steps = s% job% relax_irradiation_min_steps
      relax_to_this_irrad_flux = s% job% relax_to_this_irrad_flux
      set_to_this_irrad_flux = s% job% set_to_this_irrad_flux
      irrad_col_depth = s% job% irrad_col_depth
      relax_irradiation_max_yrs_dt = s% job% relax_irradiation_max_yrs_dt
      relax_mass_change = s% job% relax_mass_change
      relax_initial_mass_change = s% job% relax_initial_mass_change
      relax_mass_change_min_steps = s% job% relax_mass_change_min_steps
      relax_mass_change_max_yrs_dt = s% job% relax_mass_change_max_yrs_dt
      relax_mass_change_init_mdot = s% job% relax_mass_change_init_mdot
      relax_mass_change_final_mdot = s% job% relax_mass_change_final_mdot
      change_RTI_flag = s% job% change_RTI_flag
      change_initial_RTI_flag = s% job% change_initial_RTI_flag
      new_RTI_flag = s% job% new_RTI_flag
      change_RSP_flag = s% job% change_RSP_flag
      change_initial_RSP_flag = s% job% change_initial_RSP_flag
      new_RSP_flag = s% job% new_RSP_flag
      change_RSP2_flag = s% job% change_RSP2_flag
      change_initial_RSP2_flag = s% job% change_initial_RSP2_flag
      change_RSP2_flag_at_model_number = s% job% change_RSP2_flag_at_model_number
      new_RSP2_flag = s% job% new_RSP2_flag
      create_RSP2_model = s% job% create_RSP2_model
      change_w_div_wc_flag = s% job% change_w_div_wc_flag
      change_initial_w_div_wc_flag = s% job% change_initial_w_div_wc_flag
      new_w_div_wc_flag = s% job% new_w_div_wc_flag
      change_j_rot_flag = s% job% change_j_rot_flag
      change_initial_j_rot_flag = s% job% change_initial_j_rot_flag
      new_j_rot_flag = s% job% new_j_rot_flag
      
      create_RSP_model = s% job% create_RSP_model
      change_v_flag = s% job% change_v_flag
      change_initial_v_flag = s% job% change_initial_v_flag
      new_v_flag = s% job% new_v_flag
      change_D_omega_flag = s% job% change_D_omega_flag
      change_initial_D_omega_flag = s% job% change_initial_D_omega_flag
      new_D_omega_flag = s% job% new_D_omega_flag
      change_am_nu_rot_flag = s% job% change_am_nu_rot_flag
      change_initial_am_nu_rot_flag = s% job% change_initial_am_nu_rot_flag
      new_am_nu_rot_flag = s% job% new_am_nu_rot_flag
      use_D_omega_for_am_nu_rot = s% job% use_D_omega_for_am_nu_rot
      
      change_u_flag = s% job% change_u_flag
      change_initial_u_flag = s% job% change_initial_u_flag
      new_u_flag = s% job% new_u_flag
      
      change_reconstruction_flag = s% job% change_reconstruction_flag
      change_initial_reconstruction_flag = s% job% change_initial_reconstruction_flag
      new_reconstruction_flag = s% job% new_reconstruction_flag
      
      center_ye_limit_for_v_flag = s% job% center_ye_limit_for_v_flag
      change_rotation_flag = s% job% change_rotation_flag
      change_initial_rotation_flag = s% job% change_initial_rotation_flag
      new_rotation_flag = s% job% new_rotation_flag
      use_w_div_wc_flag_with_rotation = s% job% use_w_div_wc_flag_with_rotation
      use_j_rot_flag_with_rotation = s% job% use_j_rot_flag_with_rotation
      set_omega = s% job% set_omega
      set_initial_omega = s% job% set_initial_omega
      set_omega_step_limit = s% job% set_omega_step_limit
      set_near_zams_omega_steps = s% job% set_near_zams_omega_steps
      new_omega = s% job% new_omega
      set_omega_div_omega_crit = s% job% set_omega_div_omega_crit
      set_initial_omega_div_omega_crit = s% job% set_initial_omega_div_omega_crit
      set_omega_div_omega_crit_step_limit = s% job% set_omega_div_omega_crit_step_limit
      set_near_zams_omega_div_omega_crit_steps = s% job% set_near_zams_omega_div_omega_crit_steps
      new_omega_div_omega_crit = s% job% new_omega_div_omega_crit
      set_surface_rotation_v = s% job% set_surface_rotation_v
      set_initial_surface_rotation_v = s% job% set_initial_surface_rotation_v
      set_surf_rotation_v_step_limit = s% job% set_surf_rotation_v_step_limit
      set_near_zams_surface_rotation_v_steps = s% job% set_near_zams_surface_rotation_v_steps
      new_surface_rotation_v = s% job% new_surface_rotation_v
      relax_omega = s% job% relax_omega
      relax_initial_omega = s% job% relax_initial_omega
      near_zams_relax_omega = s% job% near_zams_relax_omega
      relax_omega_div_omega_crit = s% job% relax_omega_div_omega_crit
      relax_initial_omega_div_omega_crit = s% job% relax_initial_omega_div_omega_crit
      near_zams_relax_omega_div_omega_crit = s% job% near_zams_relax_omega_div_omega_crit
      relax_surface_rotation_v = s% job% relax_surface_rotation_v
      relax_initial_surface_rotation_v = s% job% relax_initial_surface_rotation_v
      near_zams_relax_initial_surface_rotation_v = s% job% near_zams_relax_initial_surface_rotation_v
      num_steps_to_relax_rotation = s% job% num_steps_to_relax_rotation
      relax_omega_max_yrs_dt = s% job% relax_omega_max_yrs_dt
      set_uniform_initial_composition = s% job% set_uniform_initial_composition
      initial_h1 = s% job% initial_h1
      initial_h2 = s% job% initial_h2
      initial_he3 = s% job% initial_he3
      initial_he4 = s% job% initial_he4
      initial_zfracs = s% job% initial_zfracs
      dump_missing_metals_into_heaviest = s% job% dump_missing_metals_into_heaviest
      relax_initial_composition = s% job% relax_initial_composition
      relax_initial_to_xaccrete = s% job% relax_initial_to_xaccrete
      relax_composition_filename = s% job% relax_composition_filename
      num_steps_to_relax_composition = s% job% num_steps_to_relax_composition
      timescale_for_relax_composition = s% job% timescale_for_relax_composition
      relax_initial_angular_momentum = s% job% relax_initial_angular_momentum
      max_steps_to_relax_angular_momentum = s% job% max_steps_to_relax_angular_momentum
      timescale_for_relax_angular_momentum = s% job% timescale_for_relax_angular_momentum
      max_dt_for_relax_angular_momentum = s% job% max_dt_for_relax_angular_momentum
      num_timescales_for_relax_angular_momentum = s% job% num_timescales_for_relax_angular_momentum
      relax_angular_momentum_filename = s% job% relax_angular_momentum_filename
      relax_angular_momentum_constant_omega_center = s% job% relax_angular_momentum_constant_omega_center
      relax_initial_entropy = s% job% relax_initial_entropy
      max_steps_to_relax_entropy = s% job% max_steps_to_relax_entropy
      timescale_for_relax_entropy = s% job% timescale_for_relax_entropy
      max_dt_for_relax_entropy = s% job% max_dt_for_relax_entropy
      num_timescales_for_relax_entropy = s% job% num_timescales_for_relax_entropy
      relax_entropy_filename = s% job% relax_entropy_filename
      get_entropy_for_relax_from_eos = s% job% get_entropy_for_relax_from_eos
      report_cell_for_xm = s% job% report_cell_for_xm
      set_to_xa_for_accretion = s% job% set_to_xa_for_accretion
      set_initial_to_xa_for_accretion = s% job% set_initial_to_xa_for_accretion
      set_nzlo = s% job% set_nzlo
      set_nzhi = s% job% set_nzhi
      change_Y = s% job% change_Y
      change_initial_Y = s% job% change_initial_Y
      relax_Y = s% job% relax_Y
      relax_initial_Y = s% job% relax_initial_Y
      relax_Y_minq = s% job% relax_Y_minq
      relax_Y_maxq = s% job% relax_Y_maxq
      new_Y = s% job% new_Y
      change_Z = s% job% change_Z
      change_initial_Z = s% job% change_initial_Z
      relax_Z = s% job% relax_Z
      relax_initial_Z = s% job% relax_initial_Z
      relax_Z_minq = s% job% relax_Z_minq
      relax_Z_maxq = s% job% relax_Z_maxq
      new_Z = s% job% new_Z
      steps_to_take_before_terminate = s% job% steps_to_take_before_terminate
      stop_if_this_file_exists = s% job% stop_if_this_file_exists
      set_initial_age = s% job% set_initial_age
      initial_age = s% job% initial_age
      set_initial_model_number = s% job% set_initial_model_number
      initial_model_number = s% job% initial_model_number
      set_initial_number_retries = s% job% set_initial_number_retries
      initial_number_retries = s% job% initial_number_retries
      set_initial_dt = s% job% set_initial_dt
      limit_initial_dt = s% job% limit_initial_dt
      years_for_initial_dt = s% job% years_for_initial_dt
      seconds_for_initial_dt = s% job% seconds_for_initial_dt
      
      set_initial_cumulative_energy_error = s% job% set_initial_cumulative_energy_error
      set_cumulative_energy_error = s% job% set_cumulative_energy_error
      set_cumulative_energy_error_at_step = s% job% set_cumulative_energy_error_at_step
      set_cumulative_energy_error_each_step_if_age_less_than = s% job% set_cumulative_energy_error_each_step_if_age_less_than
      new_cumulative_energy_error = s% job% new_cumulative_energy_error
      set_cumulative_energy_error_each_relax = s% job% set_cumulative_energy_error_each_relax
      
      change_net = s% job% change_net
      change_initial_net = s% job% change_initial_net
      new_net_name = s% job% new_net_name
      change_small_net = s% job% change_small_net
      change_initial_small_net = s% job% change_initial_small_net
      new_small_net_name = s% job% new_small_net_name
      
      h_he_net = s% job% h_he_net
      co_net = s% job% co_net
      adv_net = s% job% adv_net
      adjust_abundances_for_new_isos = s% job% adjust_abundances_for_new_isos
      set_uniform_xa_from_file = s% job% set_uniform_xa_from_file
      set_uniform_initial_xa_from_file = s% job% set_uniform_initial_xa_from_file
      file_for_uniform_xa = s% job% file_for_uniform_xa
      
      mix_section = s% job% mix_section
      mix_initial_section = s% job% mix_initial_section
      mix_section_nzlo = s% job% mix_section_nzlo
      mix_section_nzhi = s% job% mix_section_nzhi
      
      T9_weaklib_full_off = s% job% T9_weaklib_full_off
      T9_weaklib_full_on = s% job% T9_weaklib_full_on
      weaklib_blend_hi_Z = s% job% weaklib_blend_hi_Z
      T9_weaklib_full_off_hi_Z = s% job% T9_weaklib_full_off_hi_Z
      T9_weaklib_full_on_hi_Z = s% job% T9_weaklib_full_on_hi_Z
      
      use_suzuki_weak_rates = s% job% use_suzuki_weak_rates
      use_3a_fl87 = s% job% use_3a_fl87
      
      use_special_weak_rates = s% job% use_special_weak_rates
      special_weak_states_file = s% job% special_weak_states_file
      special_weak_transitions_file = s% job% special_weak_transitions_file
      ion_coulomb_corrections = s% job% ion_coulomb_corrections
      electron_coulomb_corrections = s% job% electron_coulomb_corrections
      
      mix_envelope_down_to_T = s% job% mix_envelope_down_to_T
      mix_initial_envelope_down_to_T = s% job% mix_initial_envelope_down_to_T
      auto_extend_net = s% job% auto_extend_net
      
      enable_adaptive_network = s% job% enable_adaptive_network
      min_x_for_keep = s% job% min_x_for_keep
      min_x_for_n = s% job% min_x_for_n
      min_x_for_add = s% job% min_x_for_add
      max_Z_for_add = s% job% max_Z_for_add
      max_N_for_add = s% job% max_N_for_add
      max_A_for_add = s% job% max_A_for_add
      
      save_model_number = s% job% save_model_number
      save_model_filename = s% job% save_model_filename
      save_model_when_terminate = s% job% save_model_when_terminate
      required_termination_code_string = s% job% required_termination_code_string
      profile_starting_model = s% job% profile_starting_model
      profile_model_number = s% job% profile_model_number
      report_retries = s% job% report_retries
      
      net_reaction_filename = s% job% net_reaction_filename
      jina_reaclib_filename = s% job% jina_reaclib_filename
      jina_reaclib_min_T9 = s% job% jina_reaclib_min_T9
      rate_tables_dir = s% job% rate_tables_dir
      rate_cache_suffix = s% job% rate_cache_suffix
      read_extra_star_job_inlist1 = s% job% read_extra_star_job_inlist1
      extra_star_job_inlist1_name = s% job% extra_star_job_inlist1_name
      read_extra_star_job_inlist2 = s% job% read_extra_star_job_inlist2
      extra_star_job_inlist2_name = s% job% extra_star_job_inlist2_name
      read_extra_star_job_inlist3 = s% job% read_extra_star_job_inlist3
      extra_star_job_inlist3_name = s% job% extra_star_job_inlist3_name
      read_extra_star_job_inlist4 = s% job% read_extra_star_job_inlist4
      extra_star_job_inlist4_name = s% job% extra_star_job_inlist4_name
      read_extra_star_job_inlist5 = s% job% read_extra_star_job_inlist5
      extra_star_job_inlist5_name = s% job% extra_star_job_inlist5_name
      set_abundance_nzlo = s% job% set_abundance_nzlo
      set_abundance_nzhi = s% job% set_abundance_nzhi
      set_abundance = s% job% set_abundance
      set_initial_abundance = s% job% set_initial_abundance
      chem_name = s% job% chem_name
      new_frac = s% job% new_frac
      set_abundance_nzlo = s% job% set_abundance_nzlo
      set_abundance_nzhi = s% job% set_abundance_nzhi
      replace_element = s% job% replace_element
      replace_initial_element = s% job% replace_initial_element
      chem_name1 = s% job% chem_name1
      chem_name2 = s% job% chem_name2
      replace_element_nzlo = s% job% replace_element_nzlo
      replace_element_nzhi = s% job% replace_element_nzhi
      do_special_test = s% job% do_special_test
      
      save_pulse_data_for_model_number = s% job% save_pulse_data_for_model_number
      save_pulse_data_when_terminate = s% job% save_pulse_data_when_terminate
      save_pulse_data_filename = s% job% save_pulse_data_filename
      
      chem_isotopes_filename = s% job% chem_isotopes_filename
      extras_lipar = s% job% extras_lipar
      extras_lrpar = s% job% extras_lrpar
      extras_lcpar = s% job% extras_lcpar
      extras_llpar = s% job% extras_llpar
      extras_ipar = s% job% extras_ipar
      extras_rpar = s% job% extras_rpar
      extras_cpar = s% job% extras_cpar
      extras_lpar = s% job% extras_lpar
      num_special_rate_factors = s% job% num_special_rate_factors
      special_rate_factor = s% job% special_rate_factor
      filename_of_special_rate = s% job% filename_of_special_rate
      
      reaction_for_special_factor = s% job% reaction_for_special_factor
      color_num_files = s% job% color_num_files
      color_file_names = s% job% color_file_names
      color_num_colors = s% job% color_num_colors
      
      warn_run_star_extras = s% job% warn_run_star_extras
      report_garbage_collection = s% job% report_garbage_collection
      num_steps_for_garbage_collection = s% job% num_steps_for_garbage_collection
   
   end subroutine set_star_job_controls_for_writing
   
   
   subroutine do_write_star_job(s, filename, ierr)
      type (star_info), pointer :: s
      character(*), intent(in) :: filename
      integer, intent(out) :: ierr
      integer :: io
      ierr = 0
      call set_star_job_controls_for_writing(s, ierr)
      if (ierr /= 0) then
         write(*, *) 'failed to set_star_job_controls_for_writing "' // trim(filename) // '"'
         return
      end if
      open(newunit = io, file = trim(filename), action = 'write', status = 'replace', iostat = ierr)
      if (ierr /= 0) then
         write(*, *) 'failed to open "' // trim(filename) // '"'
         return
      end if
      write(io, nml = star_job, iostat = ierr)
      if (ierr /= 0) then
         write(*, *) 'failed to write "' // trim(filename) // '"'
         return
      end if
      write(*, *) 'write star_job namelist values to "' // trim(filename) // '"'
      close(io)
   end subroutine do_write_star_job
   
   
   subroutine get_star_job(s, name, val, ierr)
      use utils_lib, only : StrUpCase
      type (star_info), pointer :: s
      character(len = *), intent(in) :: name
      character(len = *), intent(out) :: val
      integer, intent(out) :: ierr
      
      character(len(name) + 1) :: upper_name
      character(len = 512) :: str
      integer :: iounit, iostat, ind, i
      
      ierr = 0
      
      ! First save current controls
      call set_star_job_controls_for_writing(s, ierr)
      if(ierr/=0) return
      
      ! Write namelist to temporay file
      open(newunit = iounit, status = 'scratch')
      write(iounit, nml = star_job)
      rewind(iounit)
      
      ! Namelists get written in captials
      upper_name = trim(StrUpCase(name)) // '='
      val = ''
      ! Search for name inside namelist
      do
         read(iounit, '(A)', iostat = iostat) str
         ind = index(trim(str), trim(upper_name))
         if(ind /= 0) then
            val = str(ind + len_trim(upper_name):len_trim(str) - 1) ! Remove final comma and starting =
            do i = 1, len(val)
               if(val(i:i)=='"') val(i:i) = ' '
            end do
            exit
         end if
         if(is_iostat_end(iostat)) exit
      end do
      
      if(len_trim(val) == 0 .and. ind==0) ierr = -1
      
      close(iounit)
   
   end subroutine get_star_job
   
   subroutine set_star_job(s, name, val, ierr)
      type (star_info), pointer :: s
      character(len = *), intent(in) :: name, val
      character(len = len(name) + len(val) + 12) :: tmp
      integer, intent(out) :: ierr
      
      ierr = 0
      
      ! First save current star_job
      call set_star_job_controls_for_writing(s, ierr)
      if(ierr/=0) return
      
      tmp = ''
      tmp = '&star_job ' // trim(name) // '=' // trim(val) // '/'
      
      ! Load into namelist
      read(tmp, nml = star_job)
      
      ! Add to star
      call store_star_job_controls(s, ierr)
      if(ierr/=0) return
   
   end subroutine set_star_job


end module star_job_ctrls_io

