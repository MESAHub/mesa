! ***********************************************************************
!
!   Copyright (C) 2010-2019  Bill Paxton & The MESA Team
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

      module star_private_def

      use star_def
      use math_lib
      use utils_lib, only: is_bad_num, is_nan

      implicit none


      real(dp), parameter :: del_cntr_rho = 1d0
      real(dp), parameter :: min_cntr_rho = 3d0
      real(dp), parameter :: no_he_ignition_limit = 0.75d0
      real(dp), parameter :: no_cntr_T_drops_limit = 6.5d0

      real(dp), parameter :: center_h_gone = 1d-3
      real(dp), parameter :: center_h_going = 1d0/3d0
      real(dp), parameter :: center_he_going = 5d-2

      real(dp), parameter :: min_Eturb = 1d-20, min_Eturb_pow_1_pt_5 = 1d-30


      contains


      subroutine star_private_def_init
         use num_def
         integer :: i, im1
         logical :: okay

         include 'formats'

         okay = .true.
         
         mlt_partial_str(1:num_mlt_partials) = ''
         mlt_partial_str(mlt_dlnd00) = 'mlt_dlnd00'
         mlt_partial_str(mlt_dlnT00) = 'mlt_dlnT00'
         mlt_partial_str(mlt_dlndm1) = 'mlt_dlndm1'
         mlt_partial_str(mlt_dlnTm1) = 'mlt_dlnTm1'
         mlt_partial_str(mlt_dlnR) = 'mlt_dlnR'
         mlt_partial_str(mlt_dL) = 'mlt_dL'
         mlt_partial_str(mlt_cv_var) = 'mlt_cv_var'
         mlt_partial_str(mlt_w_div_wc_var) = 'mlt_w_div_wc_var'

         do i=1,num_mlt_partials
            if (len_trim(mlt_partial_str(i)) == 0) then
               if (i > 1) then
                  im1 = i-1
                  write(*,2) 'missing mlt_partial_str following ' // &
                     trim(mlt_partial_str(im1)), i
               else
                  write(*,2) 'missing mlt_partial_str 1'
               end if
               okay = .false.
            end if
         end do

         termination_code_str(1:num_termination_codes) = ''

         termination_code_str(t_max_age) = 'max_age'
         termination_code_str(t_max_omega_div_omega_crit) = 'max_omega_div_omega_crit'
         termination_code_str(t_peak_burn_vconv_div_cs_limit) = 'peak_burn_vconv_div_cs_limit'
         termination_code_str(t_max_model_number) = 'max_model_number'
         termination_code_str(t_eta_center_limit) = 'eta_center_limit'
         termination_code_str(t_log_center_temp_limit) = 'log_center_temp_limit'
         termination_code_str(t_log_center_temp_lower_limit) = 'log_center_temp_lower_limit'
         termination_code_str(t_center_entropy_limit) = 'center_entropy_limit'
         termination_code_str(t_center_entropy_lower_limit) = 'center_entropy_lower_limit'
         termination_code_str(t_max_entropy_limit) = 'max_entropy_limit'
         termination_code_str(t_max_entropy_lower_limit) = 'max_entropy_lower_limit'
         termination_code_str(t_log_center_density_limit) = 'log_center_density_limit'
         termination_code_str(t_log_center_density_lower_limit) = 'log_center_density_lower_limit'
         termination_code_str(t_gamma_center_limit) = 'gamma_center_limit'
         termination_code_str(t_log_max_temp_upper_limit) = 'log_max_temp_upper_limit'
         termination_code_str(t_log_max_temp_lower_limit) = 'log_max_temp_lower_limit'
         termination_code_str(t_HB_limit) = 'HB_limit'
         termination_code_str(t_star_mass_min_limit) = 'star_mass_min_limit'
         termination_code_str(t_star_mass_max_limit) = 'star_mass_max_limit'
         
         termination_code_str(t_star_species_mass_min_limit) = 'star_species_mass_min_limit'
         termination_code_str(t_star_species_mass_max_limit) = 'star_species_mass_max_limit'
         
         termination_code_str(t_xmstar_min_limit) = 'xmstar_min_limit'
         termination_code_str(t_xmstar_max_limit) = 'xmstar_max_limit'
         termination_code_str(t_envelope_mass_limit) = 'envelope_mass_limit'
         termination_code_str(t_envelope_fraction_left_limit) = 'envelope_fraction_left_limit'

         termination_code_str(t_he_core_mass_limit) = 'he_core_mass_limit'
         termination_code_str(t_c_core_mass_limit) = 'c_core_mass_limit'
         termination_code_str(t_o_core_mass_limit) = 'o_core_mass_limit'
         termination_code_str(t_si_core_mass_limit) = 'si_core_mass_limit'
         termination_code_str(t_fe_core_mass_limit) = 'fe_core_mass_limit'
         termination_code_str(t_neutron_rich_core_mass_limit) = 'neutron_rich_core_mass_limit'

         termination_code_str(t_he_layer_mass_lower_limit) = 'he_layer_mass_lower_limit'
         termination_code_str(t_abs_diff_lg_LH_lg_Ls_limit) = 'abs_diff_lg_LH_lg_Ls_limit'
         termination_code_str(t_Teff_lower_limit) = 'Teff_lower_limit'
         termination_code_str(t_Teff_upper_limit) = 'Teff_upper_limit'
         termination_code_str(t_delta_nu_lower_limit) = 'delta_nu_lower_limit'
         termination_code_str(t_delta_nu_upper_limit) = 'delta_nu_upper_limit'
         termination_code_str(t_delta_Pg_lower_limit) = 'delta_Pg_lower_limit'
         termination_code_str(t_delta_Pg_upper_limit) = 'delta_Pg_upper_limit'
         termination_code_str(t_shock_mass_upper_limit) = 'shock_mass_upper_limit'
         termination_code_str(t_mach1_mass_upper_limit) = 'mach1_mass_upper_limit'
         termination_code_str(t_photosphere_m_sub_M_center_limit) = 'photosphere_m_sub_M_center_limit'
         termination_code_str(t_photosphere_m_lower_limit) = 'photosphere_m_lower_limit'
         termination_code_str(t_photosphere_m_upper_limit) = 'photosphere_m_upper_limit'
         termination_code_str(t_photosphere_r_lower_limit) = 'photosphere_r_lower_limit'
         termination_code_str(t_photosphere_r_upper_limit) = 'photosphere_r_upper_limit'
         termination_code_str(t_log_Teff_lower_limit) = 'log_Teff_lower_limit'
         termination_code_str(t_log_Teff_upper_limit) = 'log_Teff_upper_limit'
         termination_code_str(t_log_Tsurf_lower_limit) = 'log_Tsurf_lower_limit'
         termination_code_str(t_log_Tsurf_upper_limit) = 'log_Tsurf_upper_limit'
         termination_code_str(t_log_Psurf_lower_limit) = 'log_Psurf_lower_limit'
         termination_code_str(t_log_Psurf_upper_limit) = 'log_Psurf_upper_limit'
         termination_code_str(t_log_Dsurf_lower_limit) = 'log_Dsurf_lower_limit'
         termination_code_str(t_log_Dsurf_upper_limit) = 'log_Dsurf_upper_limit'
         termination_code_str(t_log_L_lower_limit) = 'log_L_lower_limit'
         termination_code_str(t_log_L_upper_limit) = 'log_L_upper_limit'
         termination_code_str(t_log_g_lower_limit) = 'log_g_lower_limit'
         termination_code_str(t_log_g_upper_limit) = 'log_g_upper_limit'
         termination_code_str(t_power_nuc_burn_upper_limit) = 'power_nuc_burn_upper_limit'
         termination_code_str(t_power_h_burn_upper_limit) = 'power_h_burn_upper_limit'
         termination_code_str(t_power_he_burn_upper_limit) = 'power_he_burn_upper_limit'
         termination_code_str(t_power_z_burn_upper_limit) = 'power_z_burn_upper_limit'
         termination_code_str(t_power_nuc_burn_lower_limit) = 'power_nuc_burn_lower_limit'
         termination_code_str(t_power_h_burn_lower_limit) = 'power_h_burn_lower_limit'
         termination_code_str(t_power_he_burn_lower_limit) = 'power_he_burn_lower_limit'
         termination_code_str(t_power_z_burn_lower_limit) = 'power_z_burn_lower_limit'
         termination_code_str(t_phase_of_evolution_stop) = 'phase_of_evolution_stop'
         termination_code_str(t_center_R_lower_limit) = 'center_R_lower_limit'
         termination_code_str(t_center_Ye_lower_limit) = 'center_Ye_lower_limit'
         termination_code_str(t_fe_core_infall_limit) = 'fe_core_infall_limit'
         termination_code_str(t_non_fe_core_infall_limit) = 'non_fe_core_infall_limit'
         termination_code_str(t_non_fe_core_rebound_limit) = 'non_fe_core_rebound_limit'
         termination_code_str(t_v_div_csound_max_limit) = 'v_div_csound_max_limit'
         termination_code_str(t_v_div_csound_surf_limit) = 'v_div_csound_surf_limit'
         termination_code_str(t_Pgas_div_P_limit) = 'Pgas_div_P_limit'
         termination_code_str(t_Lnuc_div_L_lower_limit) = 'Lnuc_div_L_lower_limit'
         termination_code_str(t_Lnuc_div_L_upper_limit) = 'Lnuc_div_L_upper_limit'
         termination_code_str(t_v_surf_div_v_kh_lower_limit) = 'v_surf_div_v_kh_lower_limit'
         termination_code_str(t_v_surf_div_v_kh_upper_limit) = 'v_surf_div_v_kh_upper_limit'
         termination_code_str(t_v_surf_div_v_esc_limit) = 'v_surf_div_v_esc_limit'
         termination_code_str(t_v_surf_kms_limit) = 'v_surf_kms_limit'
         termination_code_str(t_Lnuc_div_L_zams_limit) = 'Lnuc_div_L_zams_limit'
         termination_code_str(t_phase_PreMS) = 'phase_PreMS'
         termination_code_str(t_phase_ZAMS) = 'phase_ZAMS'
         termination_code_str(t_phase_IAMS) = 'phase_IAMS'
         termination_code_str(t_phase_TAMS) = 'phase_TAMS'
         termination_code_str(t_phase_He_Burn) = 'phase_He_Burn'
         termination_code_str(t_phase_ZACHeB) = 'phase_ZACHeB' 
         termination_code_str(t_phase_TACHeB) = 'phase_TACHeB'
         termination_code_str(t_phase_TP_AGB) = 'phase_TP_AGB'
         termination_code_str(t_phase_C_Burn) = 'phase_C_Burn'
         termination_code_str(t_phase_Ne_Burn) = 'phase_Ne_Burn'
         termination_code_str(t_phase_O_Burn) = 'phase_O_Burn'
         termination_code_str(t_phase_Si_Burn) = 'phase_Si_Burn'
         termination_code_str(t_phase_WDCS) = 'phase_WDCS'
         termination_code_str(t_xa_central_lower_limit) = 'xa_central_lower_limit'
         termination_code_str(t_xa_central_upper_limit) = 'xa_central_upper_limit'
         termination_code_str(t_xa_surface_lower_limit) = 'xa_surface_lower_limit'
         termination_code_str(t_xa_surface_upper_limit) = 'xa_surface_upper_limit'
         termination_code_str(t_xa_average_lower_limit) = 'xa_average_lower_limit'
         termination_code_str(t_xa_average_upper_limit) = 'xa_average_upper_limit'
         termination_code_str(t_surface_accel_div_grav_limit) = 'surface_accel_div_grav_limit'
         termination_code_str(t_adjust_mesh_failed) = 'adjust_mesh_failed'
         termination_code_str(t_dt_is_zero) = 'dt_is_zero'
         termination_code_str(t_min_timestep_limit) = 'min_timestep_limit'
         termination_code_str(t_failed_prepare_for_new_try) = 'failed_prepare_for_new_try'
         termination_code_str(t_negative_total_angular_momentum) = 'negative_total_angular_momentum'
         termination_code_str(t_max_number_retries) = 'max_number_retries'
         termination_code_str(t_redo_limit) = 'redo_limit'
         termination_code_str(t_solve_burn) = 'solve_burn'
         termination_code_str(t_solve_hydro) = 'solve_hydro'
         termination_code_str(t_solve_mix) = 'solve_mix'
         termination_code_str(t_solve_omega_mix) = 'solve_omega_mix'
         termination_code_str(t_timestep_controller) = 'timestep_controller'
         termination_code_str(t_relax_finished_okay) = 'relax_finished_okay'
         termination_code_str(t_delta_total_energy) = 'delta total energy'
         termination_code_str(t_cumulative_extra_heating_limit) = 'cumulative extra heating'
         termination_code_str(t_max_explicit_hydro_nsteps) = 'reached max explicit hydro nsteps'
         termination_code_str(t_max_period_number) = 'reached max number of periods'
         termination_code_str(t_max_abs_rel_run_E_err) = 'exceeded max abs rel_run_E_err'

         termination_code_str(t_extras_check_model) = 'extras_check_model'
         termination_code_str(t_extras_finish_step) = 'extras_finish_step'

         termination_code_str(t_xtra1) = 'customize by setting termination_code_str(t_xtra1)'
         termination_code_str(t_xtra2) = 'customize by setting termination_code_str(t_xtra2)'
         termination_code_str(t_xtra3) = 'customize by setting termination_code_str(t_xtra3)'
         termination_code_str(t_xtra4) = 'customize by setting termination_code_str(t_xtra4)'
         termination_code_str(t_xtra5) = 'customize by setting termination_code_str(t_xtra5)'
         termination_code_str(t_xtra6) = 'customize by setting termination_code_str(t_xtra6)'
         termination_code_str(t_xtra7) = 'customize by setting termination_code_str(t_xtra7)'
         termination_code_str(t_xtra8) = 'customize by setting termination_code_str(t_xtra8)'
         termination_code_str(t_xtra9) = 'customize by setting termination_code_str(t_xtra9)'

         do i=1,num_termination_codes
            if (len_trim(termination_code_str(i)) == 0) then
               if (i > 1) then
                  im1 = i-1
                  write(*,2) 'missing termination_code_str following ' // &
                     trim(termination_code_str(im1)), i
               else
                  write(*,2) 'missing termination_code_str 1'
               end if
               okay = .false.
            end if
         end do

         dt_why_str(1:numTlim) = ''

         dt_why_str(Tlim_struc) = 'varcontrol'
         dt_why_str(Tlim_max_timestep_factor) = 'max increase'
         dt_why_str(Tlim_min_timestep_factor) = 'max decrease'
         dt_why_str(Tlim_solver) = 'solver iters'
         dt_why_str(Tlim_num_burn_steps) = 'burn steps'
         dt_why_str(Tlim_num_diff_solver_steps) = 'diff steps'
         dt_why_str(Tlim_num_diff_solver_iters) = 'diff iters'
         dt_why_str(Tlim_burn_max_num_substeps) = 'burn steps'
         dt_why_str(Tlim_burn_max_num_iters) = 'burn iters'
         dt_why_str(Tlim_max_mix_fixup) = 'mix fixup'
         dt_why_str(Tlim_dH) = 'dH'
         dt_why_str(Tlim_dHe) = 'dHe'
         dt_why_str(Tlim_dHe3) = 'dHe3'
         dt_why_str(Tlim_dX) = 'dX'
         dt_why_str(Tlim_dH_div_H) = 'dH/H'
         dt_why_str(Tlim_dHe_div_He) = 'dHe/He'
         dt_why_str(Tlim_dHe3_div_He3) = 'dHe3/He3'
         dt_why_str(Tlim_dX_div_X) = 'dX/X'
         dt_why_str(Tlim_dL_div_L) = 'dL/L'
         dt_why_str(Tlim_dlgP) = 'lgP'
         dt_why_str(Tlim_dlgRho) = 'lgRho'
         dt_why_str(Tlim_dlgT) = 'lgT'
         dt_why_str(Tlim_dlgE) = 'lgE'
         dt_why_str(Tlim_dlgR) = 'lgR'
         dt_why_str(Tlim_dlgL_nuc_cat) = 'Lnuc_cat'
         dt_why_str(Tlim_dlgL_H) = 'Lnuc_H'
         dt_why_str(Tlim_dlgL_He) = 'Lnuc_He'
         dt_why_str(Tlim_dlgL_z) = 'Lnuc_z'
         dt_why_str(Tlim_dlgL_photo) = 'Lnuc_photo'
         dt_why_str(Tlim_dlgL_nuc) = 'Lnuc'
         dt_why_str(Tlim_dvsurf_kms) = 'vsurf'
         dt_why_str(Tlim_dlgTeff) = 'lgTeff'
         dt_why_str(Tlim_dlgRho_cntr) = 'lgRho_cntr'
         dt_why_str(Tlim_dlgT_max) = 'lgT_max'
         dt_why_str(Tlim_dlgT_max_at_high_T) = 'lgT_max_hi_T'
         dt_why_str(Tlim_dlgT_cntr) = 'lgT_cntr'
         dt_why_str(Tlim_dlgP_cntr) = 'lgP_cntr'
         dt_why_str(Tlim_dlog_eps_nuc) = 'log_eps_nuc'
         dt_why_str(Tlim_lg_XH_cntr) = 'lg_XH_cntr'
         dt_why_str(Tlim_dmstar) = 'delta_mstar'
         dt_why_str(Tlim_dt_div_min_dr_div_cs) = 'min_dr_div_cs'
         dt_why_str(Tlim_lgL_phot) = 'lgL_phot'
         dt_why_str(Tlim_lgL) = 'lgL'
         dt_why_str(Tlim_max_timestep) = 'max_dt'
         dt_why_str(Tlim_timestep_hold) = 'hold'
         dt_why_str(Tlim_dX_nuc_drop) = 'dX_burn'
         dt_why_str(Tlim_dX_div_X_cntr) = 'dX_div_X_cntr'
         dt_why_str(Tlim_lg_XHe_cntr) = 'lg_XHe_cntr'
         dt_why_str(Tlim_lg_XC_cntr) = 'lg_XC_cntr'
         dt_why_str(Tlim_lg_XNe_cntr) = 'lg_XNe_cntr'
         dt_why_str(Tlim_lg_XO_cntr) = 'lg_XO_cntr'
         dt_why_str(Tlim_lg_XSi_cntr) = 'lg_XSi_cntr'
         dt_why_str(Tlim_XH_cntr) = 'XH_cntr'
         dt_why_str(Tlim_XHe_cntr) = 'XHe_cntr'
         dt_why_str(Tlim_XC_cntr) = 'XC_cntr'
         dt_why_str(Tlim_XNe_cntr) = 'XNe_cntr'
         dt_why_str(Tlim_XO_cntr) = 'XO_cntr'
         dt_why_str(Tlim_XSi_cntr) = 'XSi_cntr'
         dt_why_str(Tlim_neg_X) = 'neg_mass_frac'
         dt_why_str(Tlim_bad_Xsum) = 'bad_X_sum'
         dt_why_str(Tlim_delta_HR) = 'delta_HR'
         dt_why_str(Tlim_del_mdot) = 'delta mdot'
         dt_why_str(Tlim_adjust_J_q) = 'adjust_J_q'
         dt_why_str(Tlim_delta_Ye_highT) = 'highT del Ye'
         dt_why_str(Tlim_error_rate_energy_conservation) = 'error rate'
         dt_why_str(Tlim_avg_v_residual) = 'avg v resid'
         dt_why_str(Tlim_max_abs_v_residual) = 'max v resid'
         dt_why_str(Tlim_avg_lgE_residual) = 'avg lgE resid'
         dt_why_str(Tlim_max_abs_lgE_residual) = 'max lgE resid'
         dt_why_str(Tlim_error_in_energy_conservation) = 'rel_E_err'
         dt_why_str(Tlim_retry) = 'retry'
         dt_why_str(Tlim_binary) = 'binary'
         dt_why_str(Tlim_error_other) = 'error_other'
         dt_why_str(Tlim_other_timestep_limit) = 'other_limit'

         do i=1,numTlim
            if (len_trim(dt_why_str(i)) == 0) then
               if (i > 1) then
                  im1 = i-1
                  write(*,2) 'missing dt_why_str following ' // trim(dt_why_str(im1)), i
               else
                  write(*,2) 'missing dt_why_str 1'
               end if
               okay = .false.
            end if
         end do

         if (.not. okay) stop 'star_private_def_init'

      end subroutine star_private_def_init


      subroutine alloc_star(id, ierr)
         use rates_def, only: rates_NACRE_if_available
         integer, intent(out) :: id, ierr
         integer :: i
         type (star_info), pointer :: s

         ierr = 0
         id = -1
!$omp critical (star_handle)
         call init_star_handles()
         do i = 1, max_star_handles
            if (.not. star_handles(i)% in_use) then
               star_handles(i)% in_use = .true.
               id = i
               exit
            end if
         end do
!$omp end critical (star_handle)
         if (id == -1) then
            ierr = -1
            return
         end if
         if (star_handles(id)% id /= id) then
            ierr = -1
            return
         end if
         s => star_handles(id)

      end subroutine alloc_star
      
      
      subroutine init_star_handles()
         integer :: i
      
         if (.not. have_initialized_star_handles) then
            do i = 1, max_star_handles
               star_handles(i)% id = i
               star_handles(i)% in_use = .false.
            end do
            have_initialized_star_handles = .true.
         end if
      
      end subroutine init_star_handles
      
      
      integer function find_next_star_id()
         integer :: id
         
         id = 0
!$omp critical (star_handle_next)
         if (have_initialized_star_handles) then
            do id = 1, max_star_handles
               if (star_handles(id)% in_use .eqv. .false.) exit
            end do
         end if
!$omp end critical (star_handle_next)         
      
      find_next_star_id  = id
      end function find_next_star_id
      


      subroutine free_star(s)
         type (star_info), pointer :: s
         star_handles(s% id)% in_use = .false.
      end subroutine free_star


      subroutine stardata_init( &
            my_mesa_dir, chem_isotopes_filename, &
            net_reaction_filename, jina_reaclib_filename, &
            use_suzuki_weak_rates, &
            use_special_weak_rates, special_weak_states_file, special_weak_transitions_file, &
            reaclib_min_T9_in, &
            rate_tables_dir, rates_cache_suffix, &
            ionization_file_prefix, ionization_Z1_suffix, &
            eosDT_cache_dir, eosPT_cache_dir, eosDE_cache_dir, &
            ionization_cache_dir, kap_cache_dir, rates_cache_dir, &
            color_num_files,color_file_names,color_num_colors,&
            ierr)
         use iso_fortran_env
         use colors_lib, only : colors_init
         use kap_lib, only: kap_init
         use eos_lib, only: eos_init
         use rates_lib, only: rates_init
         use rates_def, only: reaclib_min_T9
         use net_lib, only: net_init
         use ionization_lib, only: ionization_init
         use atm_lib
         use chem_lib
         use const_lib
         use const_def, only: mesa_data_dir
         use utils_lib
         use star_history_def, only: history_column_names_init
         use star_profile_def, only: profile_column_names_init
         character (len=*), intent(in) :: &
            my_mesa_dir, chem_isotopes_filename, net_reaction_filename, &
            jina_reaclib_filename, rate_tables_dir, &
            special_weak_states_file, special_weak_transitions_file, &
            rates_cache_suffix, &
            ionization_file_prefix, ionization_Z1_suffix, &
            eosDT_cache_dir, eosPT_cache_dir, eosDE_cache_dir, &
            ionization_cache_dir, kap_cache_dir, rates_cache_dir
         integer, intent(in) :: color_num_files
         character (len=*), intent(in) :: color_file_names(:)
         integer , intent(in):: color_num_colors(:)
         real(dp), intent(in) :: reaclib_min_T9_in
         logical, intent(in) :: use_suzuki_weak_rates, use_special_weak_rates
         integer, intent(out) :: ierr

         logical, parameter :: use_cache = .true.
         character (len=strlen) :: fname
         integer :: iounit

         logical, parameter :: dbg = .false.

         include 'formats'

         ierr = 0

         rate_tables_dir_for_star = rate_tables_dir
         rates_cache_suffix_for_star = rates_cache_suffix

         if (dbg) write(*,*) 'call const_init'
         call const_init(my_mesa_dir,ierr)
         if (ierr /= 0) return

         if (dbg) write(*,*) 'call math_init'
         call math_init()

         call star_private_def_init
         call result_reason_init

         call history_column_names_init(ierr)
         if (ierr /= 0) return

         call profile_column_names_init(ierr)
         if (ierr /= 0) return

         if (dbg) write(*,*) 'call chem_init'
         call chem_init(chem_isotopes_filename, ierr)
         if (ierr /= 0) return

         if (dbg) write(*,*) 'call colors_init'
         call colors_init(color_num_files,color_file_names,color_num_colors,ierr)
         if (ierr /= 0) return

         if (dbg) write(*,*) 'call eos_init'
         !write(*,*) 'eos_file_prefix "' // trim(eos_file_prefix) // '"'
         !write(*,*) 'eosDT_cache_dir "' // trim(eosDT_cache_dir) // '"'
         !write(*,*) 'eosPT_cache_dir "' // trim(eosPT_cache_dir) // '"'
         !write(*,*) 'eosDE_cache_dir "' // trim(eosDE_cache_dir) // '"'
         call eos_init( &
            eosDT_cache_dir, eosPT_cache_dir, eosDE_cache_dir, &
            use_cache, ierr)
         if (ierr /= 0) return

         if (dbg) write(*,*) 'call kap_init'
         !write(*,*) 'kap_cache_dir "' // trim(kap_cache_dir) // '"'
         call kap_init(use_cache, kap_cache_dir, ierr)
         if (ierr /= 0) return

         if (dbg) write(*,*) 'call rates_init'
         call rates_init( &
            net_reaction_filename, jina_reaclib_filename, &
            rate_tables_dir_for_star, &
            use_suzuki_weak_rates, &
            use_special_weak_rates, &
            special_weak_states_file, &
            special_weak_transitions_file, &
            rates_cache_dir, &
            ierr)

         if (ierr /= 0) return

         if (reaclib_min_T9_in > 0 .and. reaclib_min_T9_in /= reaclib_min_T9) then
            reaclib_min_T9 = reaclib_min_T9_in
            write(*,*)
            write(*,*)
            write(*,*)
            write(*,*)
            write(*,1) 'change reaclib_min_T9', reaclib_min_T9
            write(*,1) 'must clear data/rates_data/cache of old reaclib rates'
            write(*,*)
            write(*,*)
            write(*,*)
            write(*,*)
            write(*,*)
         end if

         if (dbg) write(*,*) 'call net_init'
         call net_init(ierr)
         if (ierr /= 0) return

         if (dbg) write(*,*) 'call atm_init'
         call atm_init(use_cache, ierr)
         if (ierr /= 0) return

         if (dbg) write(*,*) 'call ionization_init'
         call ionization_init( &
            ionization_file_prefix, ionization_Z1_suffix, &
            ionization_cache_dir, use_cache, ierr)
         if (ierr /= 0) return

         version_number = ''
         fname = trim(mesa_data_dir) // '/version_number'
         open(newunit=iounit, file=trim(fname), action='read', status='old', iostat=ierr)
         if (ierr /= 0) then
            write(*, *) 'failed to open ' // trim(fname)
            return
         end if

         read(iounit, '(A)', iostat=ierr) version_number

         if (ierr /= 0) then
            close(iounit)
            return
         end if

         close(iounit)

         write(*,'(1x,a,1x,a)') 'version_number', trim(version_number)

         !here we store useful information about the compiler and SDK
         call get_compiler_version(compiler_name,compiler_version_name)
         call get_mesasdk_version(mesasdk_version_name,ierr)
         call date_and_time(date=date)
         
      end subroutine stardata_init
      
      end module star_private_def

