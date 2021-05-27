! ***********************************************************************
!
!   Copyright (C) 2020  The MESA Team
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
!
! ***********************************************************************

      module simplex_search_data

      use const_def, only: dp
      
      implicit none   
         
      logical :: just_do_first_values, trace_chi2_info, &
         simplex_just_call_my_extras_check_model, &
         simplex_using_revised_max_yr_dt
      real(dp) :: simplex_revised_max_yr_dt

      logical :: include_logg_in_chi2
      real(dp) :: logg_target, logg_sigma
      
      logical :: include_logL_in_chi2
      real(dp) :: logL_target, logL_sigma

      logical :: include_Teff_in_chi2
      real(dp) :: Teff_target, Teff_sigma

      logical :: include_FeH_in_chi2
      real(dp) :: FeH_target, FeH_sigma
         
      logical :: include_logR_in_chi2
      real(dp) :: logR_target, logR_sigma
         
      logical :: include_age_in_chi2
      real(dp) :: age_target, age_sigma
      integer :: num_smaller_steps_before_age_target
      real(dp) :: dt_for_smaller_steps_before_age_target
         
      logical :: include_surface_Z_div_X_in_chi2
      real(dp) :: surface_Z_div_X_target, surface_Z_div_X_sigma
         
      logical :: include_surface_He_in_chi2
      real(dp) :: surface_He_target, surface_He_sigma
         
      logical :: include_Rcz_in_chi2
      real(dp) :: Rcz_target, Rcz_sigma
         
      logical :: include_solar_cs_rms_in_chi2, report_solar_cs_rms
      real(dp) :: solar_cs_rms_target, solar_cs_rms_sigma

      logical :: include_my_var1_in_chi2
      real(dp) :: my_var1_target, my_var1_sigma
      character (len=32) :: my_var1_name

      logical :: include_my_var2_in_chi2
      real(dp) :: my_var2_target, my_var2_sigma
      character (len=32) :: my_var2_name

      logical :: include_my_var3_in_chi2
      real(dp) :: my_var3_target, my_var3_sigma
      character (len=32) :: my_var3_name
      
      real(dp) :: Z_div_X_solar
      
      logical :: eval_chi2_at_target_age_only
      real(dp) :: min_age_for_chi2, max_age_for_chi2

      character (len=256) :: simplex_output_filename
      integer :: simplex_itermax, &
         simplex_fcn_calls_max, simplex_seed
      real(dp) :: &
         simplex_alpha, simplex_beta, &
         simplex_gamma, simplex_delta
      logical :: &
         simplex_enforce_bounds, &
         simplex_adaptive_random_search, &
         restart_simplex_from_file
      real(dp) :: &
         simplex_x_atol, &
         simplex_x_rtol, &
         simplex_chi2_tol, &
         simplex_centroid_weight_power

      logical :: Y_depends_on_Z
      real(dp) :: Y0, dYdZ

      logical :: vary_FeH, vary_Y, vary_mass, vary_alpha, vary_f_ov
      real(dp) :: first_FeH, first_Y, first_mass, first_alpha, first_f_ov
      real(dp) :: min_FeH, min_Y, min_mass, min_alpha, min_f_ov
      real(dp) :: max_FeH, max_Y, max_mass, max_alpha, max_f_ov
      
      logical :: vary_my_param1,vary_my_param2, vary_my_param3
      real(dp) :: &
         first_my_param1, first_my_param2, first_my_param3, &
         min_my_param1, min_my_param2, min_my_param3, &
         max_my_param1, max_my_param2, max_my_param3
      character (len=32) :: my_param1_name, my_param2_name, my_param3_name
      
      real(dp) :: f0_ov_div_f_ov, Lnuc_div_L_limit, chi2_limit
      
      real(dp) :: &
         max_yrs_dt_chi2_small_limit, chi2_limit_for_small_timesteps, &
         max_yrs_dt_chi2_smaller_limit, chi2_limit_for_smaller_timesteps, &
         max_yrs_dt_chi2_smallest_limit, chi2_limit_for_smallest_timesteps, &
         chi2_search_limit1, chi2_search_limit2, chi2_relative_increase_limit, &
         avg_age_sigma_limit, avg_model_number_sigma_limit
         
      integer :: min_num_samples_for_avg, max_num_samples_for_avg, &
         limit_num_chi2_too_big
      
      real(dp) :: min_age_limit, &
         sigmas_coeff_for_logg_limit, &
         sigmas_coeff_for_logL_limit, &
         sigmas_coeff_for_Teff_limit, &
         sigmas_coeff_for_logR_limit, &
         sigmas_coeff_for_surface_Z_div_X_limit, &
         sigmas_coeff_for_surface_He_limit, &
         sigmas_coeff_for_solar_Rcz_limit, &
         sigmas_coeff_for_solar_cs_rms_limit, &
         sigmas_coeff_for_my_var1_limit, &
         sigmas_coeff_for_my_var2_limit, &
         sigmas_coeff_for_my_var3_limit
            
      ! output controls
      logical :: write_best_model_data_for_each_sample
      integer :: num_digits
      character (len=256) :: sample_results_dir, &
         sample_results_prefix, sample_results_postfix
      
      integer :: model_num_digits

      logical :: write_profile_for_best_model
      character (len=256) :: best_model_profile_filename
      
      logical :: save_model_for_best_model
      character (len=256) :: best_model_save_model_filename
      
      logical :: save_info_for_last_model
      character (len=256) :: last_model_save_info_filename
      
      ! miscellaneous
      logical :: trace_limits, save_controls
      character (len=256) :: save_controls_filename
      
      real(dp) :: Y_frac_he3
      
      logical :: read_extra_simplex_search_inlist1
      character (len=256) :: extra_simplex_search_inlist1_name 
      
      logical :: read_extra_simplex_search_inlist2
      character (len=256) :: extra_simplex_search_inlist2_name 
      
      logical :: read_extra_simplex_search_inlist3
      character (len=256) :: extra_simplex_search_inlist3_name 
      
      logical :: read_extra_simplex_search_inlist4
      character (len=256) :: extra_simplex_search_inlist4_name 
      
      logical :: read_extra_simplex_search_inlist5
      character (len=256) :: extra_simplex_search_inlist5_name 

      namelist /simplex_search_controls/ &
         just_do_first_values, &
         simplex_just_call_my_extras_check_model, &
         simplex_using_revised_max_yr_dt, &
         simplex_revised_max_yr_dt, &
         
         trace_chi2_info, &         
         include_logg_in_chi2, &
         logg_target, logg_sigma, &
         
         include_logL_in_chi2, &
         logL_target, logL_sigma, &
         
         include_Teff_in_chi2, &
         Teff_target, Teff_sigma, &
         
         include_FeH_in_chi2, &
         FeH_target, FeH_sigma, &
         
         include_logR_in_chi2, &
         logR_target, logR_sigma, &
         
         include_age_in_chi2, &
         age_target, age_sigma, &
         num_smaller_steps_before_age_target, &
         dt_for_smaller_steps_before_age_target, &
         
         include_surface_Z_div_X_in_chi2, &
         surface_Z_div_X_target, surface_Z_div_X_sigma, &
         
         include_surface_He_in_chi2, &
         surface_He_target, surface_He_sigma, &
         
         include_Rcz_in_chi2, &
         Rcz_target, Rcz_sigma, &
         
         include_solar_cs_rms_in_chi2, &
         solar_cs_rms_target, solar_cs_rms_sigma, &
         report_solar_cs_rms, &
         
         include_my_var1_in_chi2, &
         my_var1_target, my_var1_sigma, my_var1_name, &
         
         include_my_var2_in_chi2, &
         my_var2_target, my_var2_sigma, my_var2_name, &
         
         include_my_var3_in_chi2, &
         my_var3_target, my_var3_sigma, my_var3_name, &
         
         Z_div_X_solar, &
         
         eval_chi2_at_target_age_only, &
         min_age_for_chi2, &
         max_age_for_chi2, &
         
         simplex_output_filename, &
         simplex_itermax, &
         simplex_fcn_calls_max, simplex_seed, &
         simplex_alpha, simplex_beta, &
         simplex_gamma, simplex_delta, &
         simplex_enforce_bounds, &
         simplex_adaptive_random_search, &
         restart_simplex_from_file, &
         simplex_x_atol, &
         simplex_x_rtol, &
         simplex_chi2_tol, &
         simplex_centroid_weight_power, &

         Y_depends_on_Z, Y0, dYdZ, &
         vary_FeH, vary_Y, vary_mass, vary_alpha, vary_f_ov, &
         first_FeH, first_Y, first_mass, first_alpha, first_f_ov, &
         min_FeH, min_Y, min_mass, min_alpha, min_f_ov, &
         max_FeH, max_Y, max_mass, max_alpha, max_f_ov, &
         f0_ov_div_f_ov, &
         Lnuc_div_L_limit, chi2_limit, &
         max_yrs_dt_chi2_small_limit, chi2_limit_for_small_timesteps, &
         max_yrs_dt_chi2_smaller_limit, chi2_limit_for_smaller_timesteps, &
         max_yrs_dt_chi2_smallest_limit, chi2_limit_for_smallest_timesteps, &
         chi2_search_limit1, chi2_search_limit2, &
         limit_num_chi2_too_big, chi2_relative_increase_limit, &
         avg_age_sigma_limit, avg_model_number_sigma_limit, &
         min_num_samples_for_avg, max_num_samples_for_avg, &
         min_age_limit, &
         sigmas_coeff_for_logg_limit, &
         sigmas_coeff_for_logL_limit, &
         sigmas_coeff_for_Teff_limit, &

         vary_my_param1,vary_my_param2, vary_my_param3, &
         first_my_param1, first_my_param2, first_my_param3, &
         min_my_param1, min_my_param2, min_my_param3, &
         max_my_param1, max_my_param2, max_my_param3, &
         my_param1_name, my_param2_name, my_param3_name, &
         
         sigmas_coeff_for_logR_limit, &
         sigmas_coeff_for_surface_Z_div_X_limit, &
         sigmas_coeff_for_surface_He_limit, &
         sigmas_coeff_for_solar_Rcz_limit, &
         sigmas_coeff_for_solar_cs_rms_limit, &
         sigmas_coeff_for_my_var1_limit, &
         sigmas_coeff_for_my_var2_limit, &
         sigmas_coeff_for_my_var3_limit, &

         write_best_model_data_for_each_sample, model_num_digits, &
         num_digits, sample_results_dir, &
         sample_results_prefix, sample_results_postfix, &
         
         write_profile_for_best_model, best_model_profile_filename, &
         save_model_for_best_model, best_model_save_model_filename, &
         save_info_for_last_model, last_model_save_info_filename, &
         trace_limits, save_controls, save_controls_filename, &
         Y_frac_he3, &

         read_extra_simplex_search_inlist1, &
         extra_simplex_search_inlist1_name, &
         read_extra_simplex_search_inlist2, &
         extra_simplex_search_inlist2_name, &
         read_extra_simplex_search_inlist3, &
         extra_simplex_search_inlist3_name, &
         read_extra_simplex_search_inlist4, &
         extra_simplex_search_inlist4_name, &
         read_extra_simplex_search_inlist5, &
         extra_simplex_search_inlist5_name
            

      integer :: sample_number, nvar, num_chi2_too_big
      
      integer :: &
         i_Y, i_FeH, i_mass, i_alpha, i_f_ov, &
         i_my_param1, i_my_param2, i_my_param3
      real(dp) :: &
         final_Y, final_FeH, final_mass, final_alpha, final_f_ov, &
         final_my_param1, final_my_param2, final_my_param3
      
      real(dp) :: initial_max_years_for_timestep
      logical :: okay_to_restart

      real(dp) :: &
         next_FeH_to_try, next_Y_to_try, &
         next_initial_h1_to_try, next_initial_he3_to_try, &
         next_initial_he4_to_try, &
         next_mass_to_try, next_alpha_to_try, next_f_ov_to_try, &
         next_my_param1_to_try, next_my_param2_to_try, next_my_param3_to_try

      character (len=256) :: inlist_simplex_fname

      real(dp) :: &
         best_chi2, &
         best_init_h1, &
         best_init_he3, &
         best_init_he4, &
         best_init_Z, &
         best_age, &
         best_radius, &
         best_logL, &
         best_Teff, &
         best_logg, &
         best_FeH, &
         best_logR, &
         best_surface_Z_div_X, &
         best_surface_He, &
         best_Rcz, &
         best_solar_cs_rms, &
         best_my_var1, &
         best_my_var2, &
         best_my_var3, &
         best_my_param1, &
         best_my_param2, &
         best_my_param3
         
      integer :: best_model_number
      integer :: max_num_samples
      integer :: scan_grid_skip_number
             
      real(dp), pointer, dimension(:) :: &
         sample_chi2, &
         sample_age, &
         sample_init_Y, &
         sample_init_FeH, &
         sample_init_h1, &
         sample_init_he3, &
         sample_init_he4, &
         sample_init_Z, &
         sample_mass, &
         sample_alpha, &
         sample_f_ov, &
         sample_radius, &
         sample_logL, &
         sample_Teff, &
         sample_logg, &
         sample_FeH, &
         sample_logR, &
         sample_surface_Z_div_X, &
         sample_surface_He, &
         sample_Rcz, &
         sample_solar_cs_rms, &
         sample_my_var1, &
         sample_my_var2, &
         sample_my_var3, &
         sample_my_param1, &
         sample_my_param2, &
         sample_my_param3
         
      integer, pointer, dimension(:) :: &
         sample_index_by_chi2, &
         sample_model_number, &
         sample_op_code
         
      real(dp) :: simplex_max_dt_next
            
      real(dp) :: avg_age_top_samples, avg_age_sigma, &
         avg_model_number_top_samples, avg_model_number_sigma

      real(dp) :: chi2, &
         initial_Y, initial_FeH, initial_Z_div_X, &
         logg, FeH, logR, surface_Z_div_X, surface_He, Rcz, solar_cs_rms, &
         my_var1, my_var2, my_var3, my_param1, my_param2, my_param3

      integer :: star_id, star_model_number
      integer :: num_chi2_terms
      
      real(dp) :: &
         current_Y, &
         current_FeH, &
         current_mass, &
         current_alpha, &
         current_f_ov, &
         current_h1, &
         current_he3, &
         current_he4, &
         current_Z

      ! solar sound speed data
      logical :: have_sound_speed_data = .false.
      integer, parameter :: npts = 79
      real(dp), dimension(npts) :: data_r, data_csound, data_width
      
      logical, parameter :: scale_simplex_params = .false. ! experimental


      ! interfaces for procedure pointers
      abstract interface

         subroutine set_my_vars_interface(id, ierr)
            integer, intent(in) :: id
            integer, intent(out) :: ierr
         end subroutine set_my_vars_interface

         subroutine will_set_my_param_interface(id, i, new_value, ierr)
            use const_def, only: dp
            integer, intent(in) :: id
            integer, intent(in) :: i ! which of my_param's will be set
            real(dp), intent(in) :: new_value
            integer, intent(out) :: ierr
         end subroutine will_set_my_param_interface

         subroutine extras_controls_interface(id, ierr)
            integer, intent(in) :: id
            integer, intent(out) :: ierr
         end subroutine extras_controls_interface

         integer function extras_check_model_interface(id)
            integer, intent(in) :: id
         end function extras_check_model_interface
         
         integer function extras_finish_step_interface(id)
            integer, intent(in) :: id
         end function extras_finish_step_interface

         subroutine extras_after_evolve_interface(id, ierr)
            integer, intent(in) :: id
            integer, intent(out) :: ierr
         end subroutine extras_after_evolve_interface

      end interface

      type simplex_procs
         procedure(set_my_vars_interface), pointer, nopass :: set_my_vars
         procedure(will_set_my_param_interface), pointer, nopass :: will_set_my_param
         procedure(extras_controls_interface), pointer, nopass :: extras_controls
         procedure(extras_check_model_interface), pointer, nopass :: extras_check_model
         procedure(extras_finish_step_interface), pointer, nopass :: extras_finish_step
         procedure(extras_after_evolve_interface), pointer, nopass :: extras_after_evolve
      end type simplex_procs

      type (simplex_procs), target, save :: star_simplex_procs
         ! gfortran seems to require "save" here.  at least it did once upon a time.


      contains


      subroutine init_simplex_search_data(ierr)
         integer, intent(out) :: ierr
         ierr = 0
         star_simplex_procs% set_my_vars => null()
         star_simplex_procs% will_set_my_param => null()
         star_simplex_procs% extras_check_model => null()
         star_simplex_procs% extras_finish_step => null()
         star_simplex_procs% extras_after_evolve => null()
      end subroutine init_simplex_search_data


      end module simplex_search_data
