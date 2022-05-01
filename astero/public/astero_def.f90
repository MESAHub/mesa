! ***********************************************************************
!
!   Copyright (C) 2020  Bill Paxton and The MESA Team
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

      module astero_def
      use star_lib
      use star_def
      use const_def
      use math_lib
      use utils_lib
      use star_pgstar
      
      implicit none
      
      ! oscillation code results
      
      integer :: num_results
      integer, pointer, dimension(:) :: el, order, em
      real(dp), pointer, dimension(:) :: inertia, cyclic_freq, growth_rate
      real(dp) :: total_time_in_oscillation_code 

         
      ! interfaces for procedure pointers
      abstract interface
      
         subroutine other_proc_interface(id, ierr)
            integer, intent(in) :: id
            integer, intent(out) :: ierr
         end subroutine other_proc_interface

         subroutine other_adipls_mode_info_interface( &
               l, order, freq, inertia, x, y, aa, data, nn, iy, iaa, ispcpr, ierr)
            use const_def, only: dp
            integer, intent(in) :: l, order
            real(dp), intent(in) :: freq, inertia
            integer, intent(in) :: nn, iy, iaa, ispcpr
            real(dp), intent(in) :: x(1:nn), y(1:iy,1:nn), aa(1:iaa,1:nn), data(8)
            integer, intent(out) :: ierr
         end subroutine other_adipls_mode_info_interface
         
      end interface
      
      type astero_info
      
         procedure(other_proc_interface), pointer, nopass :: &
            other_after_get_chi2 => null()
            
         procedure(other_adipls_mode_info_interface), pointer, nopass :: &
            other_adipls_mode_info => null()
            
      end type astero_info
      
      type (astero_info), save :: astero_other_procs
      
      logical :: use_other_after_get_chi2 = .false.
      logical :: use_other_adipls_mode_info = .false.
      

      ! chi2 = chi2_seismo*chi2_seismo_fraction &
      !      + chi2_spectroscopic_and_photometric*(1 - chi2_seismo_fraction)
      logical :: normalize_chi2_spectro
      real(dp) :: chi2_seismo_fraction

      real(dp) :: chi2_seismo_delta_nu_fraction
      real(dp) :: chi2_seismo_nu_max_fraction
      real(dp) :: chi2_seismo_r_010_fraction
      real(dp) :: chi2_seismo_r_02_fraction
      
      logical :: &
         trace_chi2_seismo_delta_nu_info, &
         trace_chi2_seismo_nu_max_info, &
         trace_chi2_seismo_ratios_info, &
         trace_chi2_seismo_frequencies_info, &
         trace_chi2_spectro_info

      logical :: &
         normalize_chi2_seismo_frequencies, &
         normalize_chi2_seismo_r_010, &
         normalize_chi2_seismo_r_02

      real(dp) :: delta_nu, delta_nu_sigma
      real(dp) :: nu_max, nu_max_sigma

      logical :: include_logg_in_chi2_spectro
      real(dp) :: logg_target, logg_sigma
      
      logical :: include_logL_in_chi2_spectro
      real(dp) :: logL_target, logL_sigma

      logical :: include_Teff_in_chi2_spectro
      real(dp) :: Teff_target, Teff_sigma

      logical :: include_FeH_in_chi2_spectro
      real(dp) :: FeH_target, FeH_sigma
         
      logical :: include_logR_in_chi2_spectro
      real(dp) :: logR_target, logR_sigma
         
      logical :: include_age_in_chi2_spectro
      real(dp) :: age_target, age_sigma
      integer :: num_smaller_steps_before_age_target
      real(dp) :: dt_for_smaller_steps_before_age_target
         
      logical :: include_surface_Z_div_X_in_chi2_spectro
      real(dp) :: surface_Z_div_X_target, surface_Z_div_X_sigma
         
      logical :: include_surface_He_in_chi2_spectro
      real(dp) :: surface_He_target, surface_He_sigma
         
      logical :: include_Rcz_in_chi2_spectro
      real(dp) :: Rcz_target, Rcz_sigma

      integer, parameter :: max_constraints = 3

      logical :: include_my_var_in_chi2_spectro(max_constraints)
      real(dp) :: my_var_target(max_constraints), my_var_sigma(max_constraints)
      character (len=32) :: my_var_name(max_constraints)
      
      real(dp) :: Z_div_X_solar

      integer, parameter :: & ! increase these if necessary
         max_nl = 1000

      ! observed modes to match to model
      integer  :: nl(0:3)
      real(dp) :: freq_target(0:3,max_nl)
      real(dp) :: freq_sigma(0:3,max_nl)
            
      character (len=100) :: search_type
      
      logical :: eval_chi2_at_target_age_only
      real(dp) :: min_age_for_chi2, max_age_for_chi2

      character (len=256) :: newuoa_output_filename
      real(dp) :: newuoa_rhoend ! search control for newuoa

      character (len=256) :: bobyqa_output_filename
      real(dp) :: bobyqa_rhoend ! search control for bobyqa
      
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
      
      character (len=256) :: scan_grid_output_filename
      logical :: restart_scan_grid_from_file
      character (len=256) :: filename_for_parameters
      integer :: max_num_from_file
      integer :: &
         file_column_for_FeH, file_column_for_Y, file_column_for_f_ov, &
         file_column_for_alpha, file_column_for_mass, &
         file_column_for_my_param1, file_column_for_my_param2, &
         file_column_for_my_param3
      character (len=256) :: from_file_output_filename
      
      logical :: Y_depends_on_Z
      real(dp) :: Y0, dYdZ

      logical :: vary_FeH, vary_Y, vary_mass, vary_alpha, vary_f_ov
      real(dp) :: first_FeH, first_Y, first_mass, first_alpha, first_f_ov
      real(dp) :: min_FeH, min_Y, min_mass, min_alpha, min_f_ov
      real(dp) :: max_FeH, max_Y, max_mass, max_alpha, max_f_ov
      real(dp) :: delta_Y, delta_FeH, delta_mass, delta_alpha, delta_f_ov
      
      logical :: vary_my_param1, vary_my_param2, vary_my_param3
      real(dp) :: &
         first_my_param1, first_my_param2, first_my_param3, &
         min_my_param1, min_my_param2, min_my_param3, &
         max_my_param1, max_my_param2, max_my_param3, &
         delta_my_param1, delta_my_param2, delta_my_param3
      character (len=32) :: my_param1_name, my_param2_name, my_param3_name
      
      real(dp) :: f0_ov_div_f_ov, Lnuc_div_L_limit, &
         chi2_spectroscopic_limit, chi2_radial_limit, chi2_delta_nu_limit
      
      real(dp) :: max_yrs_dt_when_cold, max_yrs_dt_when_warm, max_yrs_dt_when_hot, &
         max_yrs_dt_chi2_small_limit, chi2_limit_for_small_timesteps, &
         max_yrs_dt_chi2_smaller_limit, chi2_limit_for_smaller_timesteps, &
         max_yrs_dt_chi2_smallest_limit, chi2_limit_for_smallest_timesteps, &
         chi2_search_limit1, chi2_search_limit2, chi2_relative_increase_limit, &
         avg_age_sigma_limit, avg_model_number_sigma_limit
         
      integer :: min_num_samples_for_avg, max_num_samples_for_avg, &
         limit_num_chi2_too_big
      
      real(dp) :: min_age_limit
      
      real(dp) :: &
         sigmas_coeff_for_logg_limit, &
         sigmas_coeff_for_logL_limit, &
         sigmas_coeff_for_Teff_limit, &
         sigmas_coeff_for_logR_limit, &
         sigmas_coeff_for_surface_Z_div_X_limit, &
         sigmas_coeff_for_surface_He_limit, &
         sigmas_coeff_for_Rcz_limit, &
         sigmas_coeff_for_delta_nu_limit, &
         sigmas_coeff_for_my_var_limit(max_constraints)
      
      character(len=32) :: correction_scheme, &
         surf_coef1_name, surf_coef2_name
      
      real(dp) :: correction_b, correction_factor
      integer :: l0_n_obs(max_nl)
      
      ! frequency ratios for observations
      integer :: ratios_n, ratios_l0_first, ratios_l1_first
      real(dp), dimension(max_nl) :: &
         ratios_r01, sigmas_r01, &
         ratios_r10, sigmas_r10, &
         ratios_r02, sigmas_r02
      
      ! output controls
      character (len=256) :: astero_results_directory

      character (len=strlen) :: astero_results_dbl_format
      character (len=strlen) :: astero_results_int_format
      character (len=strlen) :: astero_results_txt_format

      logical :: write_best_model_data_for_each_sample
      integer :: num_digits
      character (len=256) :: sample_results_prefix, sample_results_postfix
      
      integer :: model_num_digits

      logical :: write_fgong_for_each_model
      character (len=256) :: fgong_prefix, fgong_postfix      
      logical :: write_fgong_for_best_model
      character (len=256) :: best_model_fgong_filename
      
      logical :: write_gyre_for_each_model
      character (len=256) :: gyre_prefix, gyre_postfix     
      logical :: write_gyre_for_best_model
      character (len=256) :: best_model_gyre_filename
      integer :: max_num_gyre_points
      
      logical :: write_profile_for_best_model
      character (len=256) :: best_model_profile_filename
      
      logical :: save_model_for_best_model
      character (len=256) :: best_model_save_model_filename
      
      logical :: save_info_for_last_model
      character (len=256) :: last_model_save_info_filename
      
      
      ! miscellaneous

      logical :: save_next_best_at_higher_frequency, &
         save_next_best_at_lower_frequency

      logical :: trace_limits

      logical :: save_controls
      character (len=256) :: save_controls_filename
      
      real(dp) :: Y_frac_he3
      
      integer :: save_mode_model_number = -1
      character (len=256) :: save_mode_filename
      integer :: el_to_save = -1
      integer :: order_to_save = -1
      integer :: em_to_save = -1
      
      character (len=256) :: &
         oscillation_code, &
         gyre_input_file
      logical :: gyre_non_ad
      
      logical :: trace_time_in_oscillation_code 
      
      logical :: add_atmosphere     
      logical :: keep_surface_point      
      logical :: add_center_point

      logical :: do_redistribute_mesh
      ! note: number of zones for redistribute is set in the redistrb.c input file
               
      integer :: iscan_factor(0:3) ! iscan for adipls = this factor times expected number of modes
      real(dp) :: nu_lower_factor, nu_upper_factor
         ! frequency range for adipls is set from observed frequencies times these            
      integer :: & ! misc adipls parameters
         adipls_irotkr, adipls_nprtkr, adipls_igm1kr, adipls_npgmkr

      
      logical :: read_extra_astero_search_inlist1
      character (len=256) :: extra_astero_search_inlist1_name 
      
      logical :: read_extra_astero_search_inlist2
      character (len=256) :: extra_astero_search_inlist2_name 
      
      logical :: read_extra_astero_search_inlist3
      character (len=256) :: extra_astero_search_inlist3_name 
      
      logical :: read_extra_astero_search_inlist4
      character (len=256) :: extra_astero_search_inlist4_name 
      
      logical :: read_extra_astero_search_inlist5
      character (len=256) :: extra_astero_search_inlist5_name 


      namelist /astero_search_controls/ &
         normalize_chi2_spectro, &
         chi2_seismo_fraction, &
         chi2_seismo_delta_nu_fraction, &
         chi2_seismo_nu_max_fraction, &
         chi2_seismo_r_010_fraction, &
         chi2_seismo_r_02_fraction, &
         trace_chi2_seismo_ratios_info, &
         trace_chi2_seismo_frequencies_info, &
         trace_chi2_spectro_info, &
         trace_chi2_seismo_delta_nu_info, &
         trace_chi2_seismo_nu_max_info, &
         normalize_chi2_seismo_frequencies, &
         normalize_chi2_seismo_r_010, &
         normalize_chi2_seismo_r_02, &
         delta_nu, delta_nu_sigma, &
         nu_max, nu_max_sigma, &
         
         include_logg_in_chi2_spectro, &
         logg_target, logg_sigma, &
         
         include_logL_in_chi2_spectro, &
         logL_target, logL_sigma, &
         
         include_Teff_in_chi2_spectro, &
         Teff_target, Teff_sigma, &
         
         include_FeH_in_chi2_spectro, &
         FeH_target, FeH_sigma, &
         
         include_logR_in_chi2_spectro, &
         logR_target, logR_sigma, &
         
         include_age_in_chi2_spectro, &
         age_target, age_sigma, &
         num_smaller_steps_before_age_target, &
         dt_for_smaller_steps_before_age_target, &
         
         include_surface_Z_div_X_in_chi2_spectro, &
         surface_Z_div_X_target, surface_Z_div_X_sigma, &
         
         include_surface_He_in_chi2_spectro, &
         surface_He_target, surface_He_sigma, &
         
         include_Rcz_in_chi2_spectro, &
         Rcz_target, Rcz_sigma, &

         include_my_var_in_chi2_spectro, &
         my_var_target, my_var_sigma, my_var_name, &
         
         Z_div_X_solar, &
         nl, &
         freq_target, &
         freq_sigma, &
         
         search_type, &
         
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

         newuoa_output_filename, &
         newuoa_rhoend, &
         bobyqa_output_filename, &
         bobyqa_rhoend, &
         
         scan_grid_output_filename, &
         restart_scan_grid_from_file, &
         filename_for_parameters, &
         max_num_from_file, &
         file_column_for_FeH, file_column_for_Y, file_column_for_f_ov, &
         file_column_for_alpha, file_column_for_mass, &
         file_column_for_my_param1, file_column_for_my_param2, &
         file_column_for_my_param3, &
         from_file_output_filename, &
         Y_depends_on_Z, Y0, dYdZ, &
         vary_FeH, vary_Y, vary_mass, vary_alpha, vary_f_ov, &
         vary_my_param1, vary_my_param2, vary_my_param3, &
         first_FeH, first_Y, first_mass, first_alpha, first_f_ov, &
         first_my_param1, first_my_param2, first_my_param3, &
         min_FeH, min_Y, min_mass, min_alpha, min_f_ov, &
         min_my_param1, min_my_param2, min_my_param3, &
         max_FeH, max_Y, max_mass, max_alpha, max_f_ov, &
         max_my_param1, max_my_param2, max_my_param3, &
         delta_Y, delta_FeH, delta_mass, delta_alpha, delta_f_ov, &
         delta_my_param1, delta_my_param2, delta_my_param3, &
         my_param1_name, my_param2_name, my_param3_name, &
         f0_ov_div_f_ov, &
         Lnuc_div_L_limit, chi2_spectroscopic_limit, &
         chi2_radial_limit, chi2_delta_nu_limit, &
         max_yrs_dt_when_cold, max_yrs_dt_when_warm, max_yrs_dt_when_hot, &
         max_yrs_dt_chi2_small_limit, chi2_limit_for_small_timesteps, &
         max_yrs_dt_chi2_smaller_limit, chi2_limit_for_smaller_timesteps, &
         max_yrs_dt_chi2_smallest_limit, chi2_limit_for_smallest_timesteps, &
         chi2_search_limit1, chi2_search_limit2, &
         limit_num_chi2_too_big, chi2_relative_increase_limit, &
         avg_age_sigma_limit, avg_model_number_sigma_limit, &
         min_num_samples_for_avg, max_num_samples_for_avg, &
         sigmas_coeff_for_logg_limit, &
         sigmas_coeff_for_logL_limit, &
         sigmas_coeff_for_Teff_limit, &
         
         sigmas_coeff_for_logR_limit, &
         sigmas_coeff_for_surface_Z_div_X_limit, &
         sigmas_coeff_for_surface_He_limit, &
         sigmas_coeff_for_Rcz_limit, &
         sigmas_coeff_for_my_var_limit, &
         
         sigmas_coeff_for_delta_nu_limit, &
         min_age_limit, &
         correction_scheme, surf_coef1_name, surf_coef2_name, &
         correction_b, correction_factor, &
         l0_n_obs, &

         astero_results_directory, &

         astero_results_dbl_format, &
         astero_results_int_format, &
         astero_results_txt_format, &
         write_best_model_data_for_each_sample, &
         num_digits, &
         sample_results_prefix, sample_results_postfix, &
         model_num_digits, &
         
         write_fgong_for_each_model, &
         fgong_prefix, fgong_postfix, &
         write_fgong_for_best_model, best_model_fgong_filename, &
         
         write_gyre_for_each_model, &
         gyre_prefix, gyre_postfix, &
         write_gyre_for_best_model, best_model_gyre_filename, &
         max_num_gyre_points, &
         
         write_profile_for_best_model, best_model_profile_filename, &
         save_model_for_best_model, best_model_save_model_filename, &
         save_info_for_last_model, last_model_save_info_filename, &
         trace_limits, save_controls, save_controls_filename, &
         Y_frac_he3, &
         save_mode_model_number, save_mode_filename, &
         save_next_best_at_higher_frequency, &
         save_next_best_at_lower_frequency, &
         
         oscillation_code, &
         gyre_input_file, &
         gyre_non_ad, &
         
         el_to_save, &
         order_to_save, &
         em_to_save, &
         add_atmosphere, &
         keep_surface_point, &
         add_center_point, &
         do_redistribute_mesh, &
         trace_time_in_oscillation_code, &
         iscan_factor, &
         adipls_irotkr, adipls_nprtkr, adipls_igm1kr, adipls_npgmkr, &
         nu_lower_factor, nu_upper_factor, &
         read_extra_astero_search_inlist1, &
         extra_astero_search_inlist1_name, &
         read_extra_astero_search_inlist2, &
         extra_astero_search_inlist2_name, &
         read_extra_astero_search_inlist3, &
         extra_astero_search_inlist3_name, &
         read_extra_astero_search_inlist4, &
         extra_astero_search_inlist4_name, &
         read_extra_astero_search_inlist5, &
         extra_astero_search_inlist5_name
            
      
      ! pgstar plots

      logical :: echelle_win_flag, echelle_file_flag
      integer :: echelle_file_interval
      character (len=256) :: echelle_file_dir, echelle_file_prefix, &
         echelle_best_model_file_prefix, echelle_title
      real :: &
         echelle_win_width, echelle_win_aspect_ratio, &
         echelle_file_width, echelle_file_aspect_ratio, &
         echelle_xleft, echelle_xright, echelle_ybot, echelle_ytop, &
         echelle_txt_scale, echelle_delta_nu, echelle_model_alt_y_shift
      logical :: &
         show_echelle_next_best_at_higher_frequency, &
         show_echelle_next_best_at_lower_frequency, &
         show_echelle_annotation1, &
         show_echelle_annotation2, &
         show_echelle_annotation3      

      logical :: ratios_win_flag, ratios_file_flag
      integer :: ratios_file_interval
      character (len=256) :: ratios_file_dir, ratios_file_prefix, &
         ratios_best_model_file_prefix, ratios_title
      real :: &
         ratios_win_width, ratios_win_aspect_ratio, &
         ratios_file_width, ratios_file_aspect_ratio, &
         ratios_xleft, ratios_xright, ratios_ybot, &
         ratios_ytop, ratios_txt_scale, ratios_margin_sig_factor
      logical :: &
         show_ratios_annotation1, &
         show_ratios_annotation2, &
         show_ratios_annotation3      
      
      logical :: read_extra_astero_pgstar_inlist1
      character (len=256) :: extra_astero_pgstar_inlist1_name 
      
      logical :: read_extra_astero_pgstar_inlist2
      character (len=256) :: extra_astero_pgstar_inlist2_name 
      
      logical :: read_extra_astero_pgstar_inlist3
      character (len=256) :: extra_astero_pgstar_inlist3_name 
      
      logical :: read_extra_astero_pgstar_inlist4
      character (len=256) :: extra_astero_pgstar_inlist4_name 
      
      logical :: read_extra_astero_pgstar_inlist5
      character (len=256) :: extra_astero_pgstar_inlist5_name 
         
      namelist /astero_pgstar_controls/ &
         echelle_win_flag, echelle_file_flag, &
         echelle_file_interval, &
         echelle_file_dir, echelle_file_prefix, echelle_best_model_file_prefix, &
         echelle_win_width, echelle_win_aspect_ratio, &
         echelle_file_width, echelle_file_aspect_ratio, &
         echelle_xleft, echelle_xright, echelle_ybot, echelle_ytop, &
         echelle_txt_scale, echelle_delta_nu, echelle_title, &
         echelle_model_alt_y_shift, &
         show_echelle_next_best_at_higher_frequency, &
         show_echelle_next_best_at_lower_frequency, &
         show_echelle_annotation1, &
         show_echelle_annotation2, &
         show_echelle_annotation3, &
         ratios_win_flag, ratios_file_flag, &
         ratios_file_interval, ratios_title, &
         ratios_file_dir, ratios_file_prefix, ratios_best_model_file_prefix, &
         ratios_win_width, ratios_win_aspect_ratio, &
         ratios_file_width, ratios_file_aspect_ratio, &
         ratios_xleft, ratios_xright, ratios_ybot, &
         ratios_ytop, ratios_txt_scale, ratios_margin_sig_factor, &
         show_ratios_annotation1, &
         show_ratios_annotation2, &
         show_ratios_annotation3, &
         read_extra_astero_pgstar_inlist1, extra_astero_pgstar_inlist1_name, &
         read_extra_astero_pgstar_inlist2, extra_astero_pgstar_inlist2_name, &
         read_extra_astero_pgstar_inlist3, extra_astero_pgstar_inlist3_name, &
         read_extra_astero_pgstar_inlist4, extra_astero_pgstar_inlist4_name, &
         read_extra_astero_pgstar_inlist5, extra_astero_pgstar_inlist5_name



      ! private data
      
      
      ! working storage for models and search results
      real(dp) :: model_freq(0:3,max_nl)
      real(dp) :: model_freq_corr(0:3,max_nl)
      real(dp) :: model_inertia(0:3,max_nl)
      integer  :: model_order(0:3,max_nl)

      ! there are no alt frequencies for l=0 but we need the data so
      ! we can pass arrays containing data for all l

      ! next best fit at higher frequency
      real(dp) :: model_freq_alt_up(0:3,max_nl)
      real(dp) :: model_freq_corr_alt_up(0:3,max_nl)
      real(dp) :: model_inertia_alt_up(0:3,max_nl)
      integer  :: model_order_alt_up(0:3,max_nl)
         
      ! next best fit at lower frequency
      real(dp) :: model_freq_alt_down(0:3,max_nl)
      real(dp) :: model_freq_corr_alt_down(0:3,max_nl)
      real(dp) :: model_inertia_alt_down(0:3,max_nl)
      integer  :: model_order_alt_down(0:3,max_nl)

      integer, dimension(max_nl) :: i2_for_r02

      ! frequency ratios for model
      integer :: model_ratios_n, model_ratios_l0_first, model_ratios_l1_first
      real(dp), dimension(max_nl) :: &
         model_ratios_r01, &
         model_ratios_r10, &
         model_ratios_r02
      
      logical :: have_radial, have_nonradial
      
      real(dp) :: min_sample_chi2_so_far = -1
      integer :: sample_number, nvar, num_chi2_too_big
      
      integer :: &
         i_Y, i_FeH, &
         i_mass, i_alpha, i_f_ov, &
         i_my_param1, i_my_param2, i_my_param3
      real(dp) :: &
         final_Y, final_FeH, &
         final_mass, final_alpha, final_f_ov, &
         final_my_param1, final_my_param2, final_my_param3

      real(dp) :: initial_max_years_for_timestep
      logical :: okay_to_restart
      real(dp) :: nu_max_sun, delta_nu_sun

      real(dp) :: &
         next_FeH_to_try, next_Y_to_try, &
         next_initial_h1_to_try, next_initial_he3_to_try, &
         next_initial_he4_to_try, &
         next_mass_to_try, next_alpha_to_try, next_f_ov_to_try, &
         next_my_param1_to_try, next_my_param2_to_try, next_my_param3_to_try

      real(dp) :: avg_nu_obs, avg_radial_n
      real(dp) :: chi2_seismo_freq_fraction

      character (len=256) :: inlist_astero_fname

      real(dp) :: &
         best_chi2, &
         best_chi2_seismo, &
         best_chi2_spectro, &
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
         best_my_param1, &
         best_my_param2, &
         best_my_param3, &
         best_delta_nu, &
         best_nu_max, &
         best_surf_coef1, &
         best_surf_coef2

      real(dp) :: best_my_var(max_constraints)
         
      integer :: &
         best_model_number
         
      integer  :: best_order(0:3,max_nl)
      real(dp) :: best_freq(0:3,max_nl)
      real(dp) :: best_freq_corr(0:3,max_nl)
      real(dp) :: best_inertia(0:3,max_nl)

      real(dp), dimension(max_nl) :: &
         best_ratios_r01, &
         best_ratios_r10, &
         best_ratios_r02
         
      integer :: max_num_samples
      integer :: scan_grid_skip_number
             
      real(dp), pointer, dimension(:) :: &
         sample_chi2, &
         sample_chi2_seismo, &
         sample_chi2_spectro, &
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
         sample_my_param1, &
         sample_my_param2, &
         sample_my_param3, &
         sample_delta_nu, &
         sample_nu_max, &
         sample_surf_coef1, &
         sample_surf_coef2

      real(dp), pointer, dimension(:,:) :: sample_my_var
         
      integer, pointer, dimension(:) :: &
         sample_index_by_chi2, &
         sample_model_number, &
         sample_op_code
         
      integer,  pointer, dimension(:,:,:) :: sample_order
      real(dp), pointer, dimension(:,:,:) :: sample_freq
      real(dp), pointer, dimension(:,:,:) :: sample_freq_corr
      real(dp), pointer, dimension(:,:,:) :: sample_inertia

      real(dp), pointer, dimension(:,:) :: &
         sample_ratios_r01, &
         sample_ratios_r10, &
         sample_ratios_r02

      real(dp) :: astero_max_dt_next
            
      real(dp) :: avg_age_top_samples, avg_age_sigma, &
         avg_model_number_top_samples, avg_model_number_sigma

      real(dp) :: a_div_r, correction_a, correction_r, &
         a1, a3, power_law_a, power_law_b, sonoi_a, sonoi_b, &
         surf_coef1, surf_coef2, &
         nu_max_model, delta_nu_model, avg_nu_model, chi2, &
         chi2_seismo, chi2_spectro, chi2_radial, chi2_delta_nu, chi2_nu_max, &
         chi2_r_010_ratios, chi2_r_02_ratios, chi2_frequencies, &
         initial_Y, initial_FeH, initial_Z_div_X, &
         logg, FeH, logR, surface_Z_div_X, surface_He, Rcz, &
         my_var(max_constraints), my_param1, my_param2, my_param3

      integer :: star_id, star_model_number
      integer :: num_chi2_seismo_terms, num_chi2_spectro_terms
      
      ! current values for parameters set by adipls_extras_controls
      real(dp) :: &
         current_Y, &
         current_FeH, &
         current_mass, &
         current_alpha, &
         current_f_ov, &
         current_h1, &
         current_he3, &
         current_he4, &
         current_Z, &
         current_my_param1, &
         current_my_param2, &
         current_my_param3

      integer, parameter :: num_extra_history_columns = 5

      type (pgstar_win_file_data), pointer :: p_echelle, p_ratios

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

      end interface

      type astero_procs
         procedure(set_my_vars_interface), pointer, nopass :: set_my_vars
         procedure(will_set_my_param_interface), pointer, nopass :: will_set_my_param
         procedure(extras_startup_interface), pointer, nopass :: extras_startup
         procedure(extras_controls_interface), pointer, nopass :: extras_controls
         procedure(extras_check_model_interface), pointer, nopass :: extras_check_model
         procedure(extras_finish_step_interface), pointer, nopass :: extras_finish_step
         procedure(extras_after_evolve_interface), pointer, nopass :: extras_after_evolve
         procedure (how_many_extra_history_columns_interface), pointer, nopass :: &
            how_many_extra_history_columns
         procedure (data_for_extra_history_columns_interface), pointer, nopass :: &
            data_for_extra_history_columns
         procedure (how_many_extra_profile_columns_interface), pointer, nopass :: &
            how_many_extra_profile_columns
         procedure (data_for_extra_profile_columns_interface), pointer, nopass :: &
            data_for_extra_profile_columns
      end type astero_procs

      type (astero_procs), target, save :: star_astero_procs
         ! gfortran seems to require "save" here.  at least it did once upon a time.
      
      contains
      
      subroutine init_astero_def
         star_astero_procs% set_my_vars => null()
         star_astero_procs% will_set_my_param => null()
         star_astero_procs% extras_startup => null()
         star_astero_procs% extras_controls => null()
         star_astero_procs% extras_check_model => null()
         star_astero_procs% extras_finish_step => null()
         star_astero_procs% extras_after_evolve => null()
         star_astero_procs% how_many_extra_history_columns => null()
         star_astero_procs% data_for_extra_history_columns => null()
         star_astero_procs% how_many_extra_profile_columns => null()
         star_astero_procs% data_for_extra_profile_columns => null()
      end subroutine init_astero_def

      
      
      subroutine store_new_oscillation_results( &
            new_el, new_order, new_em, new_inertia, new_cyclic_freq, new_growth_rate, ierr)
         integer, intent(in) :: new_el, new_order, new_em
         real(dp), intent(in) :: new_inertia, new_cyclic_freq, new_growth_rate
         integer, intent(out) :: ierr
         
         integer :: n
         
         include 'formats'
         
         ierr = 0
         n = num_results*3/2 + 50
         if (.not. associated(el)) allocate(el(n))
         if (.not. associated(order)) allocate(order(n))
         if (.not. associated(em)) allocate(em(n))
         if (.not. associated(cyclic_freq)) allocate(cyclic_freq(n))
         if (.not. associated(growth_rate)) allocate(growth_rate(n))
         if (.not. associated(inertia)) allocate(inertia(n))
         
         if (num_results >= size(el,dim=1)) then ! enlarge
            call realloc_integer(el,n,ierr)
            if (ierr /= 0) call mesa_error(__FILE__,__LINE__)
            call realloc_integer(order,n,ierr)
            if (ierr /= 0) call mesa_error(__FILE__,__LINE__)
            call realloc_integer(em,n,ierr)
            if (ierr /= 0) call mesa_error(__FILE__,__LINE__)
            call realloc_double(cyclic_freq,n,ierr)
            if (ierr /= 0) call mesa_error(__FILE__,__LINE__)
            call realloc_double(growth_rate,n,ierr)
            if (ierr /= 0) call mesa_error(__FILE__,__LINE__)
            call realloc_double(inertia,n,ierr)
            if (ierr /= 0) call mesa_error(__FILE__,__LINE__)
         end if
         num_results = num_results+1
         
         n = num_results
         el(n) = new_el
         order(n) = new_order
         cyclic_freq(n) = new_cyclic_freq
         growth_rate(n) = new_growth_rate
         inertia(n) = new_inertia
         em(n) = new_em
         
      end subroutine store_new_oscillation_results

         
      subroutine init_sample_ptrs
         nullify( &
            sample_chi2, &
            sample_chi2_seismo, &
            sample_chi2_spectro, &
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
            sample_my_param1, &
            sample_my_param2, &
            sample_my_param3, &
            sample_radius, &
            sample_logL, &
            sample_Teff, &
            sample_logg, &
            sample_FeH, &
            sample_logR, &
            sample_surface_Z_div_X, &
            sample_surface_He, &
            sample_Rcz, &
            sample_my_var, &
            sample_delta_nu, &
            sample_nu_max, &
            sample_surf_coef1, &
            sample_surf_coef2, &
            sample_op_code, &
            sample_model_number, &
            sample_index_by_chi2, &
            sample_order, &
            sample_freq, &
            sample_freq_corr, &
            sample_inertia, &
            sample_ratios_r01, &
            sample_ratios_r10, &
            sample_ratios_r02)
      end subroutine init_sample_ptrs
      
      
      subroutine alloc_sample_ptrs(ierr)
         use utils_lib
         integer, intent(out) :: ierr
         ierr = 0
         max_num_samples = 1.5*max_num_samples + 200
         
         call realloc_double(sample_chi2,max_num_samples,ierr); if (ierr /= 0) return
         call realloc_double(sample_chi2_seismo,max_num_samples,ierr); if (ierr /= 0) return
         call realloc_double(sample_chi2_spectro,max_num_samples,ierr); if (ierr /= 0) return
         call realloc_double(sample_age,max_num_samples,ierr); if (ierr /= 0) return
         call realloc_double(sample_init_Y,max_num_samples,ierr); if (ierr /= 0) return
         call realloc_double(sample_init_FeH,max_num_samples,ierr); if (ierr /= 0) return
         call realloc_double(sample_init_h1,max_num_samples,ierr); if (ierr /= 0) return
         call realloc_double(sample_init_he3,max_num_samples,ierr); if (ierr /= 0) return
         call realloc_double(sample_init_he4,max_num_samples,ierr); if (ierr /= 0) return
         call realloc_double(sample_init_Z,max_num_samples,ierr); if (ierr /= 0) return
         call realloc_double(sample_mass,max_num_samples,ierr); if (ierr /= 0) return
         call realloc_double(sample_alpha,max_num_samples,ierr); if (ierr /= 0) return
         call realloc_double(sample_f_ov,max_num_samples,ierr); if (ierr /= 0) return
         call realloc_double(sample_my_param1,max_num_samples,ierr); if (ierr /= 0) return
         call realloc_double(sample_my_param2,max_num_samples,ierr); if (ierr /= 0) return
         call realloc_double(sample_my_param3,max_num_samples,ierr); if (ierr /= 0) return
         call realloc_double(sample_radius,max_num_samples,ierr); if (ierr /= 0) return
         call realloc_double(sample_logL,max_num_samples,ierr); if (ierr /= 0) return
         call realloc_double(sample_Teff,max_num_samples,ierr); if (ierr /= 0) return
         call realloc_double(sample_logg,max_num_samples,ierr); if (ierr /= 0) return
         call realloc_double(sample_FeH,max_num_samples,ierr); if (ierr /= 0) return
         
         call realloc_double(sample_logR,max_num_samples,ierr); if (ierr /= 0) return
         call realloc_double(sample_surface_Z_div_X,max_num_samples,ierr); if (ierr /= 0) return
         call realloc_double(sample_surface_He,max_num_samples,ierr); if (ierr /= 0) return
         call realloc_double(sample_Rcz,max_num_samples,ierr); if (ierr /= 0) return
         call realloc_double2(sample_my_var,max_constraints,max_num_samples,ierr); if (ierr /= 0) return
         
         call realloc_double(sample_delta_nu,max_num_samples,ierr); if (ierr /= 0) return
         call realloc_double(sample_nu_max,max_num_samples,ierr); if (ierr /= 0) return
         call realloc_double(sample_surf_coef1,max_num_samples,ierr); if (ierr /= 0) return
         call realloc_double(sample_surf_coef2,max_num_samples,ierr); if (ierr /= 0) return

         call realloc_integer(sample_index_by_chi2,max_num_samples,ierr); if (ierr /= 0) return
            
         call realloc_integer(sample_op_code,max_num_samples,ierr); if (ierr /= 0) return
         call realloc_integer(sample_model_number,max_num_samples,ierr); if (ierr /= 0) return

         call realloc_integer2_modes(sample_order,max_nl,max_num_samples,ierr); if (ierr /= 0) return
         call realloc_double2_modes(sample_freq,max_nl,max_num_samples,ierr); if (ierr /= 0) return
         call realloc_double2_modes(sample_freq_corr,max_nl,max_num_samples,ierr); if (ierr /= 0) return
         call realloc_double2_modes(sample_inertia,max_nl,max_num_samples,ierr); if (ierr /= 0) return

         call realloc_double2(sample_ratios_r01,max_nl,max_num_samples,ierr); if (ierr /= 0) return
         call realloc_double2(sample_ratios_r10,max_nl,max_num_samples,ierr); if (ierr /= 0) return
         call realloc_double2(sample_ratios_r02,max_nl,max_num_samples,ierr); if (ierr /= 0) return

      end subroutine alloc_sample_ptrs

   
      ! for the frequency arrays sample_{order,freq,freq_corr,inertia}, the first index
      ! is 0:3, so here are some specific realloc routines for that case
      ! basically copied from utils/public/utils_lib.f
      subroutine realloc_double2_modes(ptr,new_size1,new_size2,ierr)
         real(dp), pointer :: ptr(:,:,:)
         integer, intent(in) :: new_size1,new_size2
         integer, intent(out) :: ierr
         real(dp), pointer :: new_ptr(:,:,:)
         integer :: l, i1,i2, i, j
         ierr = 0
         allocate(new_ptr(0:3,new_size1,new_size2),stat=ierr)
         if (ierr /= 0) return
         if (associated(ptr)) then
            i1 = min(new_size1,size(ptr,dim=1))
            i2 = min(new_size2,size(ptr,dim=2))
            ! ifort uses stack for array copy temp storage
            ! for large copies, this can produce seg faults
            ! doing the explicit loops seems to be safe
            !new_ptr(0:3,1:i1,1:i2) = ptr(0:3,1:i1,1:i2)
            do l=0,3
               do i=1,i1
                  do j=1,i2
                     new_ptr(l,i,j) = ptr(l,i,j)
                  end do
               end do
            end do
            deallocate(ptr)
         end if
         ptr => new_ptr
      end subroutine realloc_double2_modes


      subroutine realloc_integer2_modes(ptr,new_size1,new_size2,ierr)
         integer, pointer :: ptr(:,:,:)
         integer, intent(in) :: new_size1,new_size2
         integer, intent(out) :: ierr
         integer, pointer :: new_ptr(:,:,:)
         integer :: l, i1, i2, i, j
         ierr = 0
         allocate(new_ptr(0:3,new_size1,new_size2),stat=ierr)
         if (ierr /= 0) return
         if (associated(ptr)) then
            i1 = min(new_size1,size(ptr,dim=1))
            i2 = min(new_size2,size(ptr,dim=2))
            ! ifort uses stack for array copy temp storage
            ! for large copies, this can produce seg faults
            ! doing the explicit loops seems to be safe
            !new_ptr(0:3,1:i1,1:i2) = ptr(0:3,1:i1,1:i2)
            do l=0,3
               do i=1,i1
                  do j=1,i2
                     new_ptr(l,i,j) = ptr(l,i,j)
                  end do
               end do
            end do
            deallocate(ptr)
         end if
         ptr => new_ptr
      end subroutine realloc_integer2_modes


      subroutine read_astero_search_controls(filename, ierr)
         character (len=*), intent(in) :: filename
         integer, intent(out) :: ierr
         ! initialize controls to default values
         include 'astero_search.defaults'
         ierr = 0
         call read1_astero_search_inlist(filename, 1, ierr)
      end subroutine read_astero_search_controls
         
         
      recursive subroutine read1_astero_search_inlist(filename, level, ierr)
         character (len=*), intent(in) :: filename
         integer, intent(in) :: level  
         integer, intent(out) :: ierr
         
         logical :: read_extra1, read_extra2, read_extra3, read_extra4, read_extra5
         character (len=256) :: message, extra1, extra2, extra3, extra4, extra5
         integer :: unit
         
         if (level >= 10) then
            write(*,*) 'ERROR: too many levels of nested extra star_job inlist files'
            ierr = -1
            return
         end if
         
         ierr = 0
         unit=alloc_iounit(ierr)
         if (ierr /= 0) return
         
         open(unit=unit, file=trim(filename), action='read', delim='quote', iostat=ierr)
         if (ierr /= 0) then
            write(*, *) 'Failed to open astero search inlist file ', trim(filename)
         else
            read(unit, nml=astero_search_controls, iostat=ierr)  
            close(unit)
            if (ierr /= 0) then
               write(*, *) &
                  'Failed while trying to read astero search inlist file ', trim(filename)
               write(*, '(a)') trim(message)
               write(*, '(a)') &
                  'The following runtime error message might help you find the problem'
               write(*, *) 
               open(unit=unit, file=trim(filename), &
                  action='read', delim='quote', status='old', iostat=ierr)
               read(unit, nml=astero_search_controls)
               close(unit)
            end if  
         end if
         call free_iounit(unit)
         if (ierr /= 0) return
         
         ! recursive calls to read other inlists
         
         read_extra1 = read_extra_astero_search_inlist1
         read_extra_astero_search_inlist1 = .false.
         extra1 = extra_astero_search_inlist1_name
         extra_astero_search_inlist1_name = 'undefined'
         
         read_extra2 = read_extra_astero_search_inlist2
         read_extra_astero_search_inlist2 = .false.
         extra2 = extra_astero_search_inlist2_name
         extra_astero_search_inlist2_name = 'undefined'
         
         read_extra3 = read_extra_astero_search_inlist3
         read_extra_astero_search_inlist3 = .false.
         extra3 = extra_astero_search_inlist3_name
         extra_astero_search_inlist3_name = 'undefined'
         
         read_extra4 = read_extra_astero_search_inlist4
         read_extra_astero_search_inlist4 = .false.
         extra4 = extra_astero_search_inlist4_name
         extra_astero_search_inlist4_name = 'undefined'
         
         read_extra5 = read_extra_astero_search_inlist5
         read_extra_astero_search_inlist5 = .false.
         extra5 = extra_astero_search_inlist5_name
         extra_astero_search_inlist5_name = 'undefined'
         
         if (read_extra1) then
            !write(*,*) 'read extra astero_search inlist1 from ' // trim(extra1)
            call read1_astero_search_inlist(extra1, level+1, ierr)
            if (ierr /= 0) return
         end if
         
         if (read_extra2) then
            !write(*,*) 'read extra astero_search inlist2 from ' // trim(extra2)
            call read1_astero_search_inlist(extra2, level+1, ierr)
            if (ierr /= 0) return
         end if
         
         if (read_extra3) then
            !write(*,*) 'read extra astero_search inlist3 from ' // trim(extra3)
            call read1_astero_search_inlist(extra3, level+1, ierr)
            if (ierr /= 0) return
         end if
         
         if (read_extra4) then
            !write(*,*) 'read extra astero_search inlist4 from ' // trim(extra4)
            call read1_astero_search_inlist(extra4, level+1, ierr)
            if (ierr /= 0) return
         end if
         
         if (read_extra5) then
            write(*,*) 'read extra astero_search inlist5 from ' // trim(extra5)
            call read1_astero_search_inlist(extra5, level+1, ierr)
            if (ierr /= 0) return
         end if
         
      end subroutine read1_astero_search_inlist


      subroutine write_astero_search_controls(filename_in, ierr)
         use utils_lib
         character(*), intent(in) :: filename_in
         integer, intent(out) :: ierr
         character (len=256) :: filename
         integer :: unit
         ierr = 0
         filename = trim(filename_in)
         if (len_trim(filename) == 0) filename = 'astero_search_controls.out'
         unit=alloc_iounit(ierr)
         if (ierr /= 0) then
            write(*,*) 'failed to alloc iounit in write_astero_search_controls'
            return
         end if
         ! NOTE: when open namelist file, must include delim='APOSTROPHE'
         open(unit=unit, file=trim(filename), action='write', delim='APOSTROPHE', iostat=ierr)
         if (ierr /= 0) then
            write(*, *) 'Failed to open ', trim(filename)
         else
            write(unit, nml=astero_search_controls)
            close(unit)
         end if
         call free_iounit(unit)
         
         write(*,'(A)')
         write(*,*) 'saved initial &astero_search_controls to ' // trim(filename)
         write(*,'(A)')
         write(*,'(A)')

      end subroutine write_astero_search_controls
   
   
      subroutine read_astero_pgstar_controls(filename, ierr)
         character (len=*), intent(in) :: filename
         integer, intent(out) :: ierr
         
         ! initialize controls to default values
         include 'astero_pgstar.defaults'
         
         ierr = 0
         call read1_astero_pgstar_inlist(filename, 1, ierr)
         
      end subroutine read_astero_pgstar_controls
      
   
      recursive subroutine read1_astero_pgstar_inlist(filename, level, ierr)
         character (len=*), intent(in) :: filename
         integer, intent(in) :: level  
         integer, intent(out) :: ierr
         
         logical :: read_extra1, read_extra2, read_extra3, read_extra4, read_extra5
         character (len=256) :: message, extra1, extra2, extra3, extra4, extra5
         integer :: unit
         
         if (level >= 10) then
            write(*,*) 'ERROR: too many levels of nested extra star_job inlist files'
            ierr = -1
            return
         end if
         
         ierr = 0
         unit=alloc_iounit(ierr)
         if (ierr /= 0) return
         
         open(unit=unit, file=trim(filename), action='read', delim='quote', iostat=ierr)
         if (ierr /= 0) then
            write(*, *) 'Failed to open astero pgstar inlist file ', trim(filename)
         else
            read(unit, nml=astero_pgstar_controls, iostat=ierr)  
            close(unit)
            if (ierr /= 0) then
               write(*, *) &
                  'Failed while trying to read astero pgstar inlist file ', trim(filename)
               write(*, '(a)') &
                  'The following runtime error message might help you find the problem'
               write(*, *) 
               open(unit=unit, file=trim(filename), &
                  action='read', delim='quote', status='old', iostat=ierr)
               read(unit, nml=astero_pgstar_controls)
               close(unit)
            end if  
         end if
         call free_iounit(unit)
         if (ierr /= 0) return
         
         ! recursive calls to read other inlists
         
         read_extra1 = read_extra_astero_pgstar_inlist1
         read_extra_astero_pgstar_inlist1 = .false.
         extra1 = extra_astero_pgstar_inlist1_name
         extra_astero_pgstar_inlist1_name = 'undefined'
         
         read_extra2 = read_extra_astero_pgstar_inlist2
         read_extra_astero_pgstar_inlist2 = .false.
         extra2 = extra_astero_pgstar_inlist2_name
         extra_astero_pgstar_inlist2_name = 'undefined'
         
         read_extra3 = read_extra_astero_pgstar_inlist3
         read_extra_astero_pgstar_inlist3 = .false.
         extra3 = extra_astero_pgstar_inlist3_name
         extra_astero_pgstar_inlist3_name = 'undefined'
         
         read_extra4 = read_extra_astero_pgstar_inlist4
         read_extra_astero_pgstar_inlist4 = .false.
         extra4 = extra_astero_pgstar_inlist4_name
         extra_astero_pgstar_inlist4_name = 'undefined'
         
         read_extra5 = read_extra_astero_pgstar_inlist5
         read_extra_astero_pgstar_inlist5 = .false.
         extra5 = extra_astero_pgstar_inlist5_name
         extra_astero_pgstar_inlist5_name = 'undefined'
         
         if (read_extra1) then
            !write(*,*) 'read extra astero_pgstar inlist1 from ' // trim(extra1)
            call read1_astero_pgstar_inlist(extra1, level+1, ierr)
            if (ierr /= 0) return
         end if
         
         if (read_extra2) then
            !write(*,*) 'read extra astero_pgstar inlist2 from ' // trim(extra2)
            call read1_astero_pgstar_inlist(extra2, level+1, ierr)
            if (ierr /= 0) return
         end if
         
         if (read_extra3) then
            !write(*,*) 'read extra astero_pgstar inlist3 from ' // trim(extra3)
            call read1_astero_pgstar_inlist(extra3, level+1, ierr)
            if (ierr /= 0) return
         end if
         
         if (read_extra4) then
            !write(*,*) 'read extra astero_pgstar inlist4 from ' // trim(extra4)
            call read1_astero_pgstar_inlist(extra4, level+1, ierr)
            if (ierr /= 0) return
         end if
         
         if (read_extra5) then
            write(*,*) 'read extra astero_pgstar inlist5 from ' // trim(extra5)
            call read1_astero_pgstar_inlist(extra5, level+1, ierr)
            if (ierr /= 0) return
         end if
         
      end subroutine read1_astero_pgstar_inlist


      subroutine save_sample_results_to_file(i_total, results_fname, ierr)
         use utils_lib
         integer, intent(in) :: i_total
         character (len=*), intent(in) :: results_fname
         integer, intent(out) :: ierr
         integer :: iounit
         write(*,*) 'save_sample_results_to_file ' // trim(results_fname)
         iounit = alloc_iounit(ierr)
         if (ierr /= 0) call mesa_error(__FILE__,__LINE__,'alloc_iounit failed')
         open(unit=iounit, file=trim(results_fname), action='write', iostat=ierr)
         if (ierr /= 0) return
         call show_all_sample_results(iounit, i_total, ierr)
         close(iounit)
         call free_iounit(iounit)         
      end subroutine save_sample_results_to_file
      
      
      subroutine set_sample_index_by_chi2
         use num_lib, only: qsort
         if (sample_number <= 0) return
         if (sample_number == 1) then
            sample_index_by_chi2(1) = 1
            return
         end if
         call qsort(sample_index_by_chi2, sample_number, sample_chi2)
      end subroutine set_sample_index_by_chi2
      
      
      subroutine show_sample_header(iounit)
         integer, intent(in) :: iounit
         
         integer :: i, j, k, l
         character (len=strlen) :: fmt
         character (len=10) :: str
         character (len=1) :: l_str

         ! column numbers
         write(fmt,'(a)') '(99' // trim(astero_results_int_format) // ')'

         k = 41 ! fixed columns

         if (chi2_seismo_fraction > 0) then

            do l = 0, 3
               k = k + 6*nl(l)
            end do

            k = k + 6*ratios_n

            if (chi2_seismo_r_02_fraction > 0) then
               k = k + 3*nl(0)
            end if

         end if

         if (search_type == 'simplex') k = k+1

         do i = 1, k
            write(iounit, fmt, advance='no') i
         end do

         write(iounit, '(a)') ! end of column numbers line

         ! column names
         write(fmt,'(a)') '(99' // trim(astero_results_txt_format) // ')'

         write(iounit, fmt, advance='no') &
            'sample', &
            'chi2', &
            'mass', &
            'init_Y', &
            'init_FeH', &
            'alpha', &
            'f_ov', &
            trim(my_param1_name), &
            trim(my_param2_name), &
            trim(my_param3_name), &
            'age', &
            
            'model_number', &
            
            'init_h1', &
            'init_he3', &
            'init_he4', &
            'init_Z', &
            'log_radius', &
            'logL', &
            'Teff', &
            'logg', &
            'Fe_H', &
            'logR', &
            'surface_Z_div_X', &
            'surface_He', &
            'Rcz', &
            trim(my_var_name(1)), &
            trim(my_var_name(2)), &
            trim(my_var_name(3)), &
            'delta_nu', &
            'nu_max', &
            trim(surf_coef1_name), &
            trim(surf_coef2_name), &
            'chi2_seismo', &
            'chi2_spectro', &
            
            'nl0', &
            'nl1', &
            'nl2', &
            'nl3', &
            'ratios_n', &
            'ratios_l0_first', &
            'ratios_l1_first'
         
         if (chi2_seismo_fraction > 0) then

            do l=0,3
               write(l_str,'(i1)') l
               do j=1,nl(l)
                  write(str,'(i3)') j
                  str = adjustl(str)
                  write(iounit, fmt, advance='no') &
                     'l' // l_str // '_order_' // trim(str), &
                     'l' // l_str // '_obs_' // trim(str), &
                     'l' // l_str // '_obs_sigma_' // trim(str), &
                     'l' // l_str // '_freq_' // trim(str), &
                     'l' // l_str // '_freq_corr_' // trim(str), &
                     'l' // l_str // '_inertia_' // trim(str)
               end do
            end do

            do j=1,ratios_n
               write(str,'(i3)') j
               str = adjustl(str)
               write(iounit, fmt, advance='no') &
                  'r01_obs_' // trim(str), &
                  'r01_obs_sigmas_' // trim(str), &
                  'r01_' // trim(str), &
                  'r10_obs_' // trim(str), &
                  'r10_obs_sigmas_' // trim(str), &
                  'r10_' // trim(str)
            end do

            if (chi2_seismo_r_02_fraction > 0) then
               do j=1,nl(0)
                  write(str,'(i3)') j
                  str = adjustl(str)
                  write(iounit, fmt, advance='no') &
                     'r02_obs_' // trim(str), &
                     'r02_obs_sigmas_' // trim(str), &
                     'r02_' // trim(str)
               end do
            end if
         
         end if

         if (search_type == 'simplex') then
            write(iounit, astero_results_txt_format, advance='no') 'step_type'
         end if
         
         write(iounit, '(a)') ! end of column names line
                                          
      end subroutine show_sample_header
      
      
      subroutine show1_sample_results(i, iounit)
         use num_lib, only: simplex_info_str
         integer, intent(in) :: i, iounit
            
         integer :: j, k, l, op_code, ierr
         character (len=256) :: info_str, fmt
         
         ierr = 0

         op_code = sample_op_code(i) 
         if (op_code <= 0) then
            info_str = ''
         else
            call simplex_info_str(op_code, info_str, ierr)
            if (ierr /= 0) then
               info_str = ''
               ierr = 0
            end if
         end if

         write(fmt,'(a)') '(' // trim(astero_results_int_format) // &
            ',10' // trim(astero_results_dbl_format) // &
            ',' // trim(astero_results_int_format) // &
            ',22' // trim(astero_results_dbl_format) // &
            ',7' // trim(astero_results_int_format) // ')'
         
         write(iounit, fmt, advance='no') i, &
            sample_chi2(i), &
            sample_mass(i), &
            sample_init_Y(i), &
            sample_init_FeH(i), &
            sample_alpha(i), &
            sample_f_ov(i), &
            sample_my_param1(i), &
            sample_my_param2(i), &
            sample_my_param3(i), &
            sample_age(i), &
            sample_model_number(i), &
            sample_init_h1(i), &
            sample_init_he3(i), &
            sample_init_he4(i), &
            sample_init_Z(i), &
            safe_log10(sample_radius(i)), &
            sample_logL(i), &
            sample_Teff(i), &
            sample_logg(i), &
            sample_FeH(i), &
            sample_logR(i), &
            sample_surface_Z_div_X(i), &
            sample_surface_He(i), &
            sample_Rcz(i), &
            sample_my_var(1,i), &
            sample_my_var(2,i), &
            sample_my_var(3,i), &
            sample_delta_nu(i), &
            sample_nu_max(i), &
            sample_surf_coef1(i), &
            sample_surf_coef2(i), &
            sample_chi2_seismo(i), &
            sample_chi2_spectro(i), &
            nl(0), &
            nl(1), &
            nl(2), &
            nl(3), &
            ratios_n, &
            ratios_l0_first, &
            ratios_l1_first
            
         if (iounit == 6) return

         if (chi2_seismo_fraction > 0) then

            write(fmt,'(a)') '(' // trim(astero_results_int_format) // &
                ',99' // trim(astero_results_dbl_format) // ')'

            do l = 0, 3
               do k = 1, nl(l)
                  write(iounit, fmt, advance='no') &
                     sample_order(l,k,i), freq_target(l,k), freq_sigma(l,k), &
                     sample_freq(l,k,i), sample_freq_corr(l,k,i), sample_inertia(l,k,i)
               end do
            end do

            write(fmt,'(a)') '(99' // trim(astero_results_dbl_format) // ')'

            do k=1,ratios_n
               write(iounit, fmt, advance='no') &
                  ratios_r01(k), sigmas_r01(k), sample_ratios_r01(k,i), &
                  ratios_r10(k), sigmas_r10(k), sample_ratios_r10(k,i)
            end do

            if (chi2_seismo_r_02_fraction > 0) then
               do k=1,nl(0)
                  write(iounit, fmt, advance='no') &
                     ratios_r02(k), sigmas_r02(k), sample_ratios_r02(k,i)
               end do
            end if
         
         end if

         if (search_type == 'simplex') then
            write(iounit, astero_results_txt_format, advance='no') trim(info_str)
         end if

         write(iounit, '(a)') ! end of line
      
      end subroutine show1_sample_results
      
      
      subroutine show_all_sample_results(iounit, i_total, ierr)
         integer, intent(in) :: iounit, i_total
         integer, intent(out) :: ierr
         integer :: i, j
         character (len=strlen) :: int_fmt, txt_fmt

         ierr = 0
         ! sort results by increasing sample_chi2
         call set_sample_index_by_chi2

         do j = 1, 3 ! line number
            i = 1 ! column number, incremented after each column is written

            call write_int('samples', sample_number)

            ! for scan_grid, we also write total size of grid
            if (i_total > 0) call write_int('total', i_total)

            call write_txt('version_number', version_number)
            call write_txt('compiler', compiler_name)
            call write_txt('build', compiler_version_name)
            call write_txt('MESA_SDK_version', mesasdk_version_name)
            call write_txt('math_backend',math_backend)
            call write_txt('date', date)
            call write_txt('search_type', search_type)

            write(iounit, '(a)') ! new line
         end do

         write(iounit, '(a)') ! blank line between header and sample data

         call show_sample_header(iounit)
         do j = 1, sample_number
            i = sample_index_by_chi2(j)
            call show1_sample_results(i, iounit)
         end do

         contains

         subroutine write_txt(name, val)
            character(len=*), intent(in) :: name, val

            select case (j)
            case (1)
               write(iounit, astero_results_int_format, advance='no') i
            case (2)
               write(iounit, astero_results_txt_format, advance='no') name
            case (3)
               write(iounit, astero_results_txt_format, advance='no') '"'//trim(val)//'"'
            end select

            i = i+1

         end subroutine write_txt

         subroutine write_int(name, val)
            character(len=*), intent(in) :: name
            integer, intent(in) :: val

            select case (j)
            case (1)
               write(iounit, astero_results_int_format, advance='no') i
            case (2)
               write(iounit, astero_results_txt_format, advance='no') name
            case (3)
               write(iounit, astero_results_int_format, advance='no') val
            end select

            i = i+1

         end subroutine write_int

      end subroutine show_all_sample_results
      
      
      subroutine show_best_el_info(io)
         integer, intent(in) :: io
         
         real(dp) :: chi2term
         integer :: i, l

         ! elaborate shenanigans to preserve header format
         character(len=8), dimension(5) :: header
          
         do l = 0, 3
            if (nl(l) > 0) then
               write(header(1), '(a2,i1)') 'l=', l
               write(header(2), '(a1,i1,a5)') 'l', l, '_freq'
               write(header(3), '(a1,i1,a5)') 'l', l, '_corr'
               write(header(4), '(a1,i1,a4)') 'l', l, '_obs'
               write(header(5), '(a1,i1,a6)') 'l', l, '_sigma'
               write(io,'(/,2a6,99a20)') &
                  trim(header(1)), 'n', 'chi2term', trim(header(2)), trim(header(3)), &
                  trim(header(4)), trim(header(5)), 'log E'

               do i = 1, nl(l)
                  if (freq_target(l,i) < 0) cycle
                  chi2term = pow2((best_freq_corr(l,i) - freq_target(l,i))/freq_sigma(l,i))
                  write(io,'(6x,i6,e20.10,99f20.10)') &
                     best_order(l,i), chi2term, best_freq(l,i), best_freq_corr(l,i), &
                     freq_target(l,i), freq_sigma(l,i), safe_log10(best_inertia(l,i))
               end do
            end if
         end do
      
      end subroutine show_best_el_info
      
      
      subroutine show_best_r010_ratios_info(io)
         integer, intent(in) :: io
         
         real(dp) :: chi2term
         integer :: i, l0_first, l1_first

         l0_first = ratios_l0_first
         l1_first = ratios_l1_first
          
         write(io,'(/,2a6,99a20)') &
            'r01', 'l=0 n', 'chi2term', 'r01', 'r01_obs', 'r01_sigma', 'l0_obs'
         do i=1,ratios_n
            chi2term = &
               pow2((model_ratios_r01(i) - ratios_r01(i))/sigmas_r01(i))
            write(io,'(6x,i6,e20.10,99f20.10)') model_order(0,i + l0_first), &
               chi2term, model_ratios_r01(i), ratios_r01(i), sigmas_r01(i), &
               freq_target(0,i + l0_first)
         end do
          
         write(io,'(/,2a6,99a20)') &
            'r10', 'l=1 n', 'chi2term', 'r10', 'r10_obs', 'r10_sigma', 'l1_obs'
         do i=1,ratios_n
            chi2term = &
               pow2((model_ratios_r10(i) - ratios_r10(i))/sigmas_r10(i))
            write(io,'(6x,i6,e20.10,99f20.10)') model_order(1,i + l1_first), &
               chi2term, model_ratios_r10(i), ratios_r10(i), sigmas_r10(i), &
               freq_target(1,i + l1_first)
         end do
               
      end subroutine show_best_r010_ratios_info

          
      subroutine show_best_r02_ratios_info(io)
         integer, intent(in) :: io
         
         real(dp) :: chi2term
         integer :: i
         
         write(io,'(/,2a6,99a20)') &
            'r02', 'l=0 n', 'chi2term', 'r02', 'r02_obs', 'r02_sigma', 'l0_obs'
         do i=1,nl(0)
            if (sigmas_r02(i) == 0d0) cycle
            chi2term = &
               pow2((model_ratios_r02(i) - ratios_r02(i))/sigmas_r02(i))
            write(io,'(6x,i6,e20.10,99f20.10)') model_order(0,i), &
               chi2term, model_ratios_r02(i), ratios_r02(i), sigmas_r02(i), &
               freq_target(0,i)
         end do
               
      end subroutine show_best_r02_ratios_info
      
      
      subroutine show_best(io)
         integer, intent(in) :: io
         
         real(dp) :: chi2term
         integer :: i
         include 'formats'
         
         if (chi2_seismo_fraction > 0) then
            call show_best_el_info(io)         
            if (chi2_seismo_r_010_fraction > 0) &
               call show_best_r010_ratios_info(io)        
            if (chi2_seismo_r_02_fraction > 0) &
               call show_best_r02_ratios_info(io)
         end if

         if (Teff_sigma > 0 .and. include_Teff_in_chi2_spectro) then
            chi2term = pow2((best_Teff - Teff_target)/Teff_sigma)
            write(io,'(A)')
            call write1('Teff chi2term', chi2term)
            call write1('Teff', best_Teff)
            call write1('Teff_obs', Teff_target)
            call write1('Teff_sigma', Teff_sigma)
         end if
         
         if (logL_sigma > 0 .and. include_logL_in_chi2_spectro) then
            chi2term = pow2((best_logL - logL_target)/logL_sigma)
            write(io,'(A)')
            call write1('logL chi2term', chi2term)
            call write1('logL', best_logL)
            call write1('logL_obs', logL_target)
            call write1('logL_sigma', logL_sigma)
         end if
         
         if (logg_sigma > 0 .and. include_logg_in_chi2_spectro) then
            chi2term = pow2((best_logg - logg_target)/logg_sigma)
            write(io,'(A)')
            call write1('logg chi2term', chi2term)
            call write1('logg', best_logg)
            call write1('logg_obs', logg_target)
            call write1('logg_sigma', logg_sigma)
         end if
         
         if (FeH_sigma > 0 .and. include_FeH_in_chi2_spectro) then
            chi2term = pow2((best_FeH - FeH_target)/FeH_sigma)
            write(io,'(A)')
            call write1('FeH chi2term', chi2term)
            call write1('FeH', best_FeH)
            call write1('FeH_obs', FeH_target)
            call write1('FeH_sigma', FeH_sigma)
         end if
         
         if (logR_sigma > 0 .and. include_logR_in_chi2_spectro) then
            chi2term = pow2((best_logR - logR_target)/logR_sigma)
            write(io,'(A)')
            call write1('logR chi2term', chi2term)
            call write1('logR', best_logR)
            call write1('logR_obs', logR_target)
            call write1('logR_sigma', logR_sigma)
         end if
         
         if (age_sigma > 0 .and. include_age_in_chi2_spectro) then
            chi2term = pow2((best_age - age_target)/age_sigma)
            write(io,'(A)')
            write(io,'(a40,e20.10,99f20.10)') 'age chi2term', chi2term
            write(io,'(a40,1pes20.10)') 'age', best_age
            write(io,'(a40,1pes20.10)') 'age_target', age_target
            write(io,'(a40,1pes20.10)') 'age_sigma', age_sigma
         end if
         
         if (surface_Z_div_X_sigma > 0 .and. &
               include_surface_Z_div_X_in_chi2_spectro) then
            chi2term = &
               pow2((best_surface_Z_div_X - surface_Z_div_X_target)/surface_Z_div_X_sigma)
            write(io,'(A)')
            write(io,'(a40,e20.10,99f20.10)') 'surface_Z_div_X chi2term', chi2term
            call write1('surface_Z_div_X', best_surface_Z_div_X)
            call write1('surface_Z_div_X_obs', surface_Z_div_X_target)
            call write1('surface_Z_div_X_sigma', surface_Z_div_X_sigma)
         end if
         
         if (surface_He_sigma > 0 .and. include_surface_He_in_chi2_spectro) then
            chi2term = pow2((best_surface_He - surface_He_target)/surface_He_sigma)
            write(io,'(A)')
            call write1('surface_He chi2term', chi2term)
            call write1('surface_He', best_surface_He)
            call write1('surface_He_obs', surface_He_target)
            call write1('surface_He_sigma', surface_He_sigma)
         end if
         
         if (Rcz_sigma > 0 .and. include_Rcz_in_chi2_spectro) then
            chi2term = pow2((best_Rcz - Rcz_target)/Rcz_sigma)
            write(io,'(A)')
            call write1('Rcz chi2term', chi2term)
            call write1('Rcz', best_Rcz)
            call write1('Rcz_obs', Rcz_target)
            call write1('Rcz_sigma', Rcz_sigma)
         end if

         do i = 1, max_constraints
            if (my_var_sigma(i) > 0 .and. include_my_var_in_chi2_spectro(i)) then
               chi2term = pow2( &
                     (best_my_var(i) - my_var_target(i))/my_var_sigma(i))
               write(io,'(A)')
               call write1(trim(my_var_name(i)) // ' chi2term', chi2term)
               call write1(trim(my_var_name(i)), best_my_var(i))
               call write1(trim(my_var_name(i)) // '_obs', my_var_target(i))
               call write1(trim(my_var_name(i)) // '_sigma', my_var_sigma(i))
            end if
         end do
         
         write(io,'(A)')
         call write1('R/Rsun', best_radius)
         call write1('logL/Lsun', best_logL)
         call write1('Teff', best_Teff)
         call write1('logg', best_logg)
         call write1('FeH', best_FeH)
         call write1('logR', best_logR)
         call write1('surface_Z_div_X', best_surface_Z_div_X)
         call write1('surface_He', best_surface_He)
         call write1('Rcz', best_Rcz)
         call write1('delta_nu', best_delta_nu)
         call write1('nu_max', best_nu_max)
         write(io,*)        
         write(io,'(a40,1pes20.10)') trim(surf_coef1_name), best_surf_coef1
         write(io,'(a40,1pes20.10)') trim(surf_coef2_name), best_surf_coef2
         write(io,*)        
         call write1('initial h1', current_h1)
         call write1('initial he3', current_he3)
         call write1('initial he4', current_he4)
         call write1('initial Y', current_Y)
         call write1('initial Z', current_Z)
         call write1('initial FeH', current_FeH)
         write(io,*)        
         call write1('mass/Msun', current_mass)
         call write1('alpha', current_alpha)
         call write1('f_ov', current_f_ov)
         call write1(trim(my_param1_name), current_my_param1)
         call write1(trim(my_param2_name), current_my_param2)
         call write1(trim(my_param3_name), current_my_param3)
         write(io,'(a40,1pes20.10)') 'age', best_age
         write(io,'(A)')
         if (chi2_seismo_fraction == 1d0) then
            call write1('chi^2 seismo', best_chi2_seismo)
         else if (chi2_seismo_fraction == 0d0) then
            call write1('chi^2 spectro', best_chi2_spectro)
         else
            call write1('chi^2 combined', best_chi2)
            call write1('chi2_seismo_fraction', chi2_seismo_fraction)
            call write1('chi^2 seismo', best_chi2_seismo)
            call write1('chi^2 spectro', best_chi2_spectro)
         end if
         write(io,'(A)')
         write(io,'(a40,i16)') 'model number', best_model_number
         write(io,'(A)')
         write(io,'(A)')
         
         contains
         
         subroutine write1(str,x)
            character (len=*), intent(in) :: str
            real(dp), intent(in) :: x
            if (abs(x) < 1d6) then
               write(io,'(a40,99f20.10)') trim(str), x
            else
               write(io,'(a40,99e20.10)') trim(str), x
            end if
         end subroutine write1

      end subroutine show_best


      subroutine read_samples_from_file(results_fname, ierr)
         use utils_lib
         character (len=*), intent(in) :: results_fname
         integer, intent(out) :: ierr
         integer :: iounit, num, i, j, model_number
         character (len=strlen) :: line
         
         include 'formats'
         
         ierr = 0         
         write(*,*) 'read samples from file ' // trim(results_fname)
         
         iounit = alloc_iounit(ierr)
         if (ierr /= 0) call mesa_error(__FILE__,__LINE__,'alloc_iounit failed')
         open(unit=iounit, file=trim(results_fname), action='read', status='old', iostat=ierr)
         if (ierr /= 0) then
            write(*,*) 'failed to open ' // trim(results_fname)
            call free_iounit(iounit) 
            return
         end if

         read(iounit, fmt='(a)') line
         read(iounit, fmt='(a)') line
         
         read(iounit, fmt=astero_results_int_format, iostat=ierr) num
         if (ierr /= 0) then
            write(*,*) 'failed to read number of samples on line 3 of ' // trim(results_fname)
            call done
            return
         end if
         
         write(*,2) 'number of samples in file', num

         do j = 4, 6
            read(iounit, fmt='(a)', iostat=ierr) line
            if (ierr /= 0) then
               write(*,'(a,i1,a)') 'failed to line ', j, ' of ' // trim(results_fname)
               write(*,'(a)') 'line <' // trim(line) // '>'
               call done
               return
            end if
         end do
         
         do while (max_num_samples < num)
            call alloc_sample_ptrs(ierr)
            if (ierr /= 0) then
               write(*,*) 'ERROR -- failed to allocate for samples'
               call done
               return
            end if
         end do
         
         do j = 1, num
            call read1_sample_from_file(j, iounit, ierr)
            if (ierr /= 0) then
               write(*,2) 'ERROR -- failed while reading sample on line', j+2
               call done
               return
            end if
         end do
                  
         sample_number = num
         write(*,2) 'number of samples read from file', num
         
         call done

         
         contains
         
         
         subroutine done
            close(iounit)
            call free_iounit(iounit)         
         end subroutine done
         

      end subroutine read_samples_from_file
      
      
      subroutine read1_sample_from_file(j, iounit, ierr)
         use num_lib, only: simplex_op_code
         integer, intent(in) :: j, iounit
         integer, intent(out) :: ierr
            
         integer :: i, k, l
         character (len=256) :: info_str, fmt
         real(dp) :: logR
         
         include 'formats'
         
         ierr = 0
         read(iounit,fmt=astero_results_int_format,advance='no',iostat=ierr) i
         if (ierr /= 0) return
         if (i <= 0 .or. i > size(sample_chi2,dim=1)) then
            write(*,2) 'invalid sample number', i
            ierr = -1
            return
         end if

         write(fmt,'(a)') '(10' // trim(astero_results_dbl_format) // &
            ',' // trim(astero_results_int_format) // &
            ',22' // trim(astero_results_dbl_format) // &
            ',7' // trim(astero_results_int_format) // ')'
         
         read(iounit, fmt, advance='no', iostat=ierr) &
            sample_chi2(i), &
            sample_mass(i), &
            sample_init_Y(i), &
            sample_init_FeH(i), &
            sample_alpha(i), &
            sample_f_ov(i), &
            sample_my_param1(i), &
            sample_my_param2(i), &
            sample_my_param3(i), &
            sample_age(i), &
            sample_model_number(i), &
            sample_init_h1(i), &
            sample_init_he3(i), &
            sample_init_he4(i), &
            sample_init_Z(i), &
            logR, &
            sample_logL(i), &
            sample_Teff(i), &
            sample_logg(i), &
            sample_FeH(i), &
            sample_logR(i), &
            sample_surface_Z_div_X(i), &
            sample_surface_He(i), &
            sample_Rcz(i), &
            sample_my_var(1,i), &
            sample_my_var(2,i), &
            sample_my_var(3,i), &
            sample_delta_nu(i), &
            sample_nu_max(i), &
            sample_surf_coef1(i), &
            sample_surf_coef2(i), &
            sample_chi2_seismo(i), &
            sample_chi2_spectro(i), &
            nl(0), &
            nl(1), &
            nl(2), &
            nl(3), &
            ratios_n, &
            ratios_l0_first, &
            ratios_l1_first
         if (failed('results')) return
            
         sample_radius(i) = exp10(logR)

         if (chi2_seismo_fraction > 0) then

            write(fmt,'(a)') '(' // trim(astero_results_int_format) // &
                ',99' // trim(astero_results_dbl_format) // ')'

            do l = 0, 3
               do k = 1, nl(l)
                  read(iounit, fmt, advance='no', iostat=ierr) &
                     sample_order(l,k,i), freq_target(l,k), freq_sigma(l,k), &
                     sample_freq(l,k,i), sample_freq_corr(l,k,i), sample_inertia(l,k,i)
                  if (failed('freqs')) return
               end do
            end do

            write(fmt,'(a)') '(99' // trim(astero_results_dbl_format) // ')'

            do k=1,ratios_n
               read(iounit, fmt, advance='no', iostat=ierr) &
                  ratios_r01(k), sigmas_r01(k), sample_ratios_r01(k,i), &
                  ratios_r10(k), sigmas_r10(k), sample_ratios_r10(k,i)
               if (failed('ratios_r010')) return
            end do

            if (chi2_seismo_r_02_fraction > 0.0_dp) then
               do k=1,nl(0)
                  read(iounit, fmt, advance='no', iostat=ierr) &
                     ratios_r02(k), sigmas_r02(k), sample_ratios_r02(k,i)
                  if (failed('ratios_r02')) return
               end do
            end if
         
         end if
            
         read(iounit, '(a12)', iostat=ierr) info_str
         if (ierr /= 0) then
            ierr = 0
            sample_op_code(i) = 0
            return
         end if
      
         if (len_trim(info_str) == 0) then
            sample_op_code(i) = 0
         else
            sample_op_code(i) = simplex_op_code(info_str, ierr)
            if (ierr /= 0) then
               ierr = 0
               sample_op_code(i) = 0
               return
            end if
         end if
         
         
         contains
         
         
         logical function failed(str)
            character (len=*), intent(in) :: str
            include 'formats'
            failed = .false.
            if (ierr == 0) return
            write(*,2) 'failed reading ' // trim(str) // ' data for sample number', i
            failed = .true.
         end function failed
         
      
      end subroutine read1_sample_from_file


      end module astero_def

