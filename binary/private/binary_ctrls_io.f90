! ***********************************************************************
!
!   Copyright (C) 2010  Pablo Marchant
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

      module binary_ctrls_io
      
      use const_def
      use binary_def

      implicit none
      
      include "binary_controls.inc"      
      
      logical, dimension(max_extra_inlists) :: read_extra_binary_controls_inlist
      character (len=strlen), dimension(max_extra_inlists) :: extra_binary_controls_inlist_name
      
      namelist /binary_controls/ &
         ! specifications for starting model
         m1, &
         m2, &
         initial_period_in_days, &
         initial_separation_in_Rsuns, &
         initial_eccentricity, &

         ! controls for output
         history_name, &
         history_interval, &
         append_to_star_history, &
         log_directory, &
         history_dbl_format, &
         history_int_format, &
         history_txt_format, &
         photo_interval, &
         photo_digits, &
         photo_directory, &
         terminal_interval, &
         write_header_frequency, &
         extra_binary_terminal_output_file, &

         ! timestep controls
         fm, &
         fm_hard, &
         fa, &
         fa_hard, &
         fr, &
         fr_hard, &
         fj, &
         fj_hard, &
         fe, &
         fe_hard, &
         fm_limit, &
         fr_limit, &
         fe_limit, &
         fr_dt_limit, &
         fdm, &
         fdm_hard, &
         dt_softening_factor, &
         varcontrol_case_a, &
         varcontrol_case_b, &
         varcontrol_ms, &
         varcontrol_post_ms, &
         dt_reduction_factor_for_j, &
         
         ! when to stop
         accretor_overflow_terminate, &
         terminate_if_initial_overflow, &
         terminate_if_L2_overflow, &
             
         ! mass transfer controls
         mass_transfer_alpha, &
         mass_transfer_beta, &
         mass_transfer_delta, &
         mass_transfer_gamma, &
         limit_retention_by_mdot_edd, &
         use_es_opacity_for_mdot_edd, &
         use_this_for_mdot_edd_eta, &
         use_radiation_corrected_transfer_rate, &
         initial_bh_spin, &
         use_this_for_mdot_edd, &
         mdot_scheme, &
         cur_mdot_frac, &
         max_explicit_abs_mdot, &
         max_tries_to_achieve, &
         solver_type, &
         implicit_scheme_tolerance, &
         implicit_scheme_tiny_factor, &
         initial_change_factor, &
         change_factor_fraction, &
         implicit_lambda, &
         max_change_factor, &
         min_change_factor, &
         num_tries_for_increase_change_factor, &
         change_factor_increase, &
         starting_mdot, &
         roche_min_mdot, &
         min_mdot_for_implicit, &
         max_implicit_abs_mdot, &
         report_rlo_solver_progress, &
         do_enhance_wind_1, &
         do_enhance_wind_2, &
         tout_B_wind_1, &
         tout_B_wind_2, &
         do_wind_mass_transfer_1, &
         do_wind_mass_transfer_2, &
         wind_BH_alpha_1, &
         wind_BH_alpha_2, &
         wind_BH_beta_1, &
         wind_BH_beta_2, &
         max_wind_transfer_fraction_1, &
         max_wind_transfer_fraction_2, &

         ! orbital jdot controls
         do_jdot_gr, &
         do_jdot_ml, &
         do_jdot_ls, &
         do_jdot_missing_wind, &
         do_jdot_mb, &
         include_accretor_mb, &
         magnetic_braking_gamma, &
         keep_mb_on, &
         jdot_mb_min_qconv_env, &
         jdot_mb_max_qconv_env, &
         jdot_mb_max_qconv_core, &
         jdot_mb_qlim_for_check_rad_core, &
         jdot_mb_qlim_for_check_conv_env, &
         jdot_mb_scale_for_low_qconv_env, &
         jdot_mb_mass_frac_for_scale, &
         jdot_multiplier, &

         ! rotation and sync controls
         do_j_accretion, &
         do_tidal_sync, &
         sync_type_1, &
         sync_type_2, &
         sync_mode_1, &
         sync_mode_2, &
         Ftid_1, &
         Ftid_2, &
         do_initial_orbit_sync_1, &
         do_initial_orbit_sync_2, &
         tidal_reduction, &
         
         ! eccentricity controls
         do_tidal_circ, &
         circ_type_1, &
         circ_type_2, &
         use_eccentricity_enhancement, &
         max_abs_edot_tidal, &
         max_abs_edot_enhance, &
         min_eccentricity, &
         max_eccentricity, &
         anomaly_steps, &

         ! irradiation controls
         accretion_powered_irradiation, &
         col_depth_for_eps_extra, &
         use_accretor_luminosity_for_irrad, &
         irrad_flux_at_std_distance, &
         std_distance_for_irradiation, &
         max_F_irr, &

         !common envelope controls
         CE_alpha, &
         CE_alpha_th, &
         CE_alpha_core, &
         CE_mass_loss_rate_high, &
         CE_mass_loss_rate_low, &
         CE_rel_rlo_for_detachment, &
         CE_years_detached_to_terminate, &
         CE_begin_at_max_implicit_abs_mdot, &
         CE_xa_diff_to_terminate, &
         CE_terminate_when_core_overflows, &
         CE_min_period_in_minutes, &
         CE_energy_factor_HII_toHI, &
         CE_energy_factor_HeII_toHeI, &
         CE_energy_factor_HeIII_toHeII, &
         CE_energy_factor_H2, &
         CE_fixed_lambda, &
         
         ! miscellaneous controls
         keep_donor_fixed, &
         mdot_limit_donor_switch, &
         use_other_rlo_mdot, &
         use_other_check_implicit_rlo, &
         use_other_implicit_function_to_solve, &
         use_other_tsync, &
         use_other_sync_spin_to_orbit, &
         use_other_mdot_edd, &
         use_other_adjust_mdots, &
         use_other_accreted_material_j, &
         use_other_jdot_gr, &
         use_other_jdot_ml, &
         use_other_jdot_ls, &
         use_other_jdot_missing_wind, &
         use_other_jdot_mb, &
         use_other_extra_jdot, &
         use_other_binary_wind_transfer, &
         use_other_edot_tidal, &
         use_other_edot_enhance, &
         use_other_extra_edot, &
         use_other_CE_init, &
         use_other_CE_rlo_mdot, &
         use_other_CE_binary_evolve_step, &
         use_other_CE_binary_finish_step, &

         ! extra files
         read_extra_binary_controls_inlist, extra_binary_controls_inlist_name

      contains
      
      
      subroutine do_one_binary_setup(b, inlist, ierr)
         use utils_lib
         type (binary_info), pointer :: b
         character (len=*), intent(in) :: inlist
         integer, intent(out) :: ierr

         include 'formats'

         call set_default_binary_controls
         call read_binary_controls(b, inlist, ierr)

         ! open additional file for binary output
         if (len_trim(b% extra_binary_terminal_output_file) /= 0) then
            open(newunit=b% extra_binary_terminal_iounit, file=trim(b% extra_binary_terminal_output_file), &
                     action='write', status='replace',iostat=ierr)
            if (ierr /= 0) then
               write(*,*) 'failed to open ' // trim(b% extra_binary_terminal_output_file)
               return
            end if
         end if

      end subroutine do_one_binary_setup


      subroutine read_binary_controls(b, filename, ierr)
         use utils_lib
         type (binary_info), pointer :: b
         character(*), intent(in) :: filename
         integer, intent(out) :: ierr
         
         call read_binary_controls_file(b, filename, 1, ierr)
         
      end subroutine read_binary_controls
         
         
      recursive subroutine read_binary_controls_file(b, filename, level, ierr)
         use utils_lib
         character(*), intent(in) :: filename
         type (binary_info), pointer :: b
         integer, intent(in) :: level
         integer, intent(out) :: ierr
         logical, dimension(max_extra_inlists) :: read_extra
         character (len=strlen) :: message
         character (len=strlen), dimension(max_extra_inlists) :: extra
         integer :: unit, i
         
         ierr = 0        
         
         if (level >= 10) then
            write(*,*) 'ERROR: too many levels of nested extra binary controls inlist files'
            ierr = -1
            return
         end if

         if (len_trim(filename) > 0) then
            open(newunit=unit, file=trim(filename), action='read', delim='quote', status='old', iostat=ierr)
            if (ierr /= 0) then
               write(*, *) 'Failed to open binary control namelist file ', trim(filename)
               return
            end if
            read(unit, nml=binary_controls, iostat=ierr)  
            close(unit)
            if (ierr /= 0) then
               write(*, *) 
               write(*, *) 
               write(*, *) 
               write(*, *) 
               write(*, '(a)') &
                  'Failed while trying to read binary control namelist file: ' // trim(filename)
               write(*, '(a)') &
                  'Perhaps the following runtime error message will help you find the problem.'
               write(*, *) 
               open(newunit=unit, file=trim(filename), action='read', delim='quote', status='old', iostat=ierr)
               read(unit, nml=binary_controls)
               close(unit)
               return
            end if
         end if
         
         call store_binary_controls(b, ierr)
         
         ! recursive calls to read other inlists
         do i=1, max_extra_inlists
            read_extra(i) = read_extra_binary_controls_inlist(i)
            read_extra_binary_controls_inlist(i) = .false.
            extra(i) = extra_binary_controls_inlist_name(i)
            extra_binary_controls_inlist_name(i) = 'undefined'
            
            if (read_extra(i)) then
               call read_binary_controls_file(b, extra(i), level+1, ierr)
               if (ierr /= 0) return
            end if
         end do
         
      end subroutine read_binary_controls_file


      subroutine set_default_binary_controls
         include 'binary_controls.defaults'
      end subroutine set_default_binary_controls


      subroutine store_binary_controls(b, ierr)
         use utils_lib, only: mkdir
         type (binary_info), pointer :: b
         integer, intent(out) :: ierr
         
         ierr = 0
         
         ! specifications for starting model
         b% m1 = m1
         b% m2 = m2
         b% initial_period_in_days = initial_period_in_days
         b% initial_separation_in_Rsuns = initial_separation_in_Rsuns
         b% initial_eccentricity = initial_eccentricity

         ! controls for output
         b% history_name = history_name
         b% history_interval = history_interval
         b% append_to_star_history = append_to_star_history
         b% log_directory = log_directory
         CALL mkdir(b% log_directory)
         b% history_dbl_format = history_dbl_format
         b% history_int_format = history_int_format
         b% history_txt_format = history_txt_format
         b% photo_interval = photo_interval
         b% photo_digits = photo_digits
         b% photo_directory = photo_directory
         CALL mkdir(b% photo_directory)
         b% terminal_interval = terminal_interval
         b% write_header_frequency = write_header_frequency
         b% extra_binary_terminal_output_file = extra_binary_terminal_output_file

         ! timestep controls
         b% fm = fm
         b% fm_hard = fm_hard
         b% fa = fa
         b% fa_hard = fa_hard
         b% fr = fr
         b% fr_hard = fr_hard
         b% fj = fj
         b% fj_hard = fj_hard
         b% fe = fe
         b% fe_hard = fe_hard
         b% fm_limit = fm_limit
         b% fr_limit = fr_limit
         b% fe_limit = fe_limit
         b% fr_dt_limit = fr_dt_limit
         b% fdm = fdm
         b% fdm_hard = fdm_hard
         b% dt_softening_factor = dt_softening_factor
         b% varcontrol_case_a = varcontrol_case_a
         b% varcontrol_case_b = varcontrol_case_b
         b% varcontrol_ms = varcontrol_ms
         b% varcontrol_post_ms = varcontrol_post_ms
         b% dt_reduction_factor_for_j = dt_reduction_factor_for_j

         ! when to stop
         b% accretor_overflow_terminate = accretor_overflow_terminate
         b% terminate_if_initial_overflow = terminate_if_initial_overflow
         b% terminate_if_L2_overflow = terminate_if_L2_overflow

         ! mass transfer controls
         b% mass_transfer_alpha = mass_transfer_alpha
         b% mass_transfer_beta = mass_transfer_beta
         b% mass_transfer_delta = mass_transfer_delta
         b% mass_transfer_gamma = mass_transfer_gamma
         b% limit_retention_by_mdot_edd = limit_retention_by_mdot_edd
         b% use_es_opacity_for_mdot_edd = use_es_opacity_for_mdot_edd
         b% use_this_for_mdot_edd_eta = use_this_for_mdot_edd_eta
         b% use_radiation_corrected_transfer_rate = use_radiation_corrected_transfer_rate
         b% initial_bh_spin = initial_bh_spin
         b% use_this_for_mdot_edd = use_this_for_mdot_edd
         b% mdot_scheme = mdot_scheme
         b% cur_mdot_frac = cur_mdot_frac
         b% max_explicit_abs_mdot = max_explicit_abs_mdot
         b% max_tries_to_achieve = max_tries_to_achieve
         b% solver_type = solver_type
         b% implicit_scheme_tolerance = implicit_scheme_tolerance
         b% implicit_scheme_tiny_factor = implicit_scheme_tiny_factor
         b% initial_change_factor = initial_change_factor
         b% change_factor_fraction = change_factor_fraction
         b% implicit_lambda = implicit_lambda
         b% max_change_factor = max_change_factor
         b% min_change_factor = min_change_factor
         b% num_tries_for_increase_change_factor = num_tries_for_increase_change_factor
         b% change_factor_increase = change_factor_increase
         b% starting_mdot = starting_mdot
         b% roche_min_mdot = roche_min_mdot
         b% min_mdot_for_implicit = min_mdot_for_implicit
         b% max_implicit_abs_mdot = max_implicit_abs_mdot
         b% report_rlo_solver_progress = report_rlo_solver_progress
         b% do_enhance_wind_1 = do_enhance_wind_1
         b% do_enhance_wind_2 = do_enhance_wind_2
         b% tout_B_wind_1 = tout_B_wind_1
         b% tout_B_wind_2 = tout_B_wind_2
         b% do_wind_mass_transfer_1 = do_wind_mass_transfer_1
         b% do_wind_mass_transfer_2 = do_wind_mass_transfer_2
         b% wind_BH_alpha_1 = wind_BH_alpha_1
         b% wind_BH_alpha_2 = wind_BH_alpha_2
         b% wind_BH_beta_1 = wind_BH_beta_1
         b% wind_BH_beta_2 = wind_BH_beta_2
         b% max_wind_transfer_fraction_1 = max_wind_transfer_fraction_1
         b% max_wind_transfer_fraction_2 = max_wind_transfer_fraction_2

         ! orbital jdot controls
         b% do_jdot_gr = do_jdot_gr
         b% do_jdot_ml = do_jdot_ml
         b% do_jdot_ls = do_jdot_ls
         b% do_jdot_missing_wind = do_jdot_missing_wind
         b% do_jdot_mb = do_jdot_mb
         b% include_accretor_mb = include_accretor_mb
         b% magnetic_braking_gamma = magnetic_braking_gamma
         b% keep_mb_on = keep_mb_on
         b% jdot_mb_min_qconv_env = jdot_mb_min_qconv_env
         b% jdot_mb_max_qconv_env = jdot_mb_max_qconv_env
         b% jdot_mb_max_qconv_core = jdot_mb_max_qconv_core
         b% jdot_mb_qlim_for_check_rad_core = jdot_mb_qlim_for_check_rad_core
         b% jdot_mb_qlim_for_check_conv_env = jdot_mb_qlim_for_check_conv_env
         b% jdot_mb_scale_for_low_qconv_env = jdot_mb_scale_for_low_qconv_env
         b% jdot_mb_mass_frac_for_scale = jdot_mb_mass_frac_for_scale
         b% jdot_multiplier = jdot_multiplier

         ! rotation and sync controls
         b% do_j_accretion = do_j_accretion
         b% do_tidal_sync = do_tidal_sync
         b% sync_type_1 = sync_type_1
         b% sync_type_2 = sync_type_2
         b% sync_mode_1 = sync_mode_1
         b% sync_mode_2 = sync_mode_2
         b% Ftid_1 = Ftid_1
         b% Ftid_2 = Ftid_2
         b% do_initial_orbit_sync_1 = do_initial_orbit_sync_1
         b% do_initial_orbit_sync_2 = do_initial_orbit_sync_2
         b% tidal_reduction = tidal_reduction
         
         ! eccentricity controls
         b% do_tidal_circ = do_tidal_circ
         b% circ_type_1 = circ_type_1
         b% circ_type_2 = circ_type_2
         b% use_eccentricity_enhancement = use_eccentricity_enhancement
         b% max_abs_edot_tidal = max_abs_edot_tidal
         b% max_abs_edot_enhance = max_abs_edot_enhance
         b% min_eccentricity = min_eccentricity
         b% max_eccentricity = max_eccentricity
         b% anomaly_steps = anomaly_steps

         ! irradiation controls
         b% accretion_powered_irradiation = accretion_powered_irradiation
         b% use_accretor_luminosity_for_irrad = use_accretor_luminosity_for_irrad
         b% col_depth_for_eps_extra = col_depth_for_eps_extra
         b% irrad_flux_at_std_distance = irrad_flux_at_std_distance
         b% std_distance_for_irradiation = std_distance_for_irradiation
         b% max_F_irr = max_F_irr

         !common envelope controls
         b% CE_alpha = CE_alpha
         b% CE_alpha_th = CE_alpha_th
         b% CE_alpha_core = CE_alpha_core
         b% CE_mass_loss_rate_high = CE_mass_loss_rate_high
         b% CE_mass_loss_rate_low = CE_mass_loss_rate_low
         b% CE_rel_rlo_for_detachment = CE_rel_rlo_for_detachment
         b% CE_years_detached_to_terminate = CE_years_detached_to_terminate
         b% CE_begin_at_max_implicit_abs_mdot = CE_begin_at_max_implicit_abs_mdot
         b% CE_xa_diff_to_terminate = CE_xa_diff_to_terminate
         b% CE_terminate_when_core_overflows = CE_terminate_when_core_overflows
         b% CE_min_period_in_minutes = CE_min_period_in_minutes
         b% CE_energy_factor_HII_toHI = CE_energy_factor_HII_toHI
         b% CE_energy_factor_HeII_toHeI = CE_energy_factor_HeII_toHeI
         b% CE_energy_factor_HeIII_toHeII = CE_energy_factor_HeIII_toHeII
         b% CE_energy_factor_H2 = CE_energy_factor_H2
         b% CE_fixed_lambda = CE_fixed_lambda
         
         ! miscellaneous controls
         b% keep_donor_fixed = keep_donor_fixed
         b% mdot_limit_donor_switch = mdot_limit_donor_switch
         b% use_other_rlo_mdot = use_other_rlo_mdot
         b% use_other_check_implicit_rlo = use_other_check_implicit_rlo
         b% use_other_implicit_function_to_solve = use_other_implicit_function_to_solve
         b% use_other_tsync = use_other_tsync
         b% use_other_sync_spin_to_orbit = use_other_sync_spin_to_orbit
         b% use_other_mdot_edd = use_other_mdot_edd
         b% use_other_adjust_mdots = use_other_adjust_mdots
         b% use_other_accreted_material_j = use_other_accreted_material_j
         b% use_other_jdot_gr = use_other_jdot_gr
         b% use_other_jdot_ml = use_other_jdot_ml
         b% use_other_jdot_ls = use_other_jdot_ls
         b% use_other_jdot_missing_wind = use_other_jdot_missing_wind
         b% use_other_jdot_mb = use_other_jdot_mb
         b% use_other_extra_jdot = use_other_extra_jdot
         b% use_other_binary_wind_transfer = use_other_binary_wind_transfer
         b% use_other_edot_tidal = use_other_edot_tidal
         b% use_other_edot_enhance = use_other_edot_enhance
         b% use_other_extra_edot = use_other_extra_edot
         b% use_other_CE_init = use_other_CE_init
         b% use_other_CE_rlo_mdot = use_other_CE_rlo_mdot
         b% use_other_CE_binary_evolve_step = use_other_CE_binary_evolve_step
         b% use_other_CE_binary_finish_step = use_other_CE_binary_finish_step
         
      end subroutine store_binary_controls


      subroutine set_binary_controls_for_writing(b, ierr)
         type (binary_info), pointer :: b
         integer, intent(out) :: ierr
         
         ierr = 0

         ! specifications for starting model
         m1 = b% m1
         m2 = b% m2
         initial_period_in_days = b% initial_period_in_days
         initial_separation_in_Rsuns = b% initial_separation_in_Rsuns
         initial_eccentricity = b% initial_eccentricity

         ! controls for output
         history_name = b% history_name
         history_interval = b% history_interval
         append_to_star_history = b% append_to_star_history
         log_directory = b% log_directory
         history_dbl_format = b% history_dbl_format
         history_int_format = b% history_int_format
         history_txt_format = b% history_txt_format
         photo_interval = b% photo_interval
         photo_digits = b% photo_digits
         photo_directory = b% photo_directory
         terminal_interval = b% terminal_interval
         write_header_frequency = b% write_header_frequency
         extra_binary_terminal_output_file = b% extra_binary_terminal_output_file

         ! timestep controls
         fm = b% fm
         fa = b% fa
         fr = b% fr
         fj = b% fj
         fe = b% fe
         fm_limit = b% fm_limit
         fr_limit = b% fr_limit
         fe_limit = b% fe_limit
         fr_dt_limit = b% fr_dt_limit
         fdm = b% fdm
         fdm_hard = b% fdm_hard
         dt_softening_factor = b% dt_softening_factor
         varcontrol_case_a = b% varcontrol_case_a
         varcontrol_case_b = b% varcontrol_case_b
         varcontrol_ms = b% varcontrol_ms
         varcontrol_post_ms = b% varcontrol_post_ms
         dt_reduction_factor_for_j = b% dt_reduction_factor_for_j

         ! when to stop
         accretor_overflow_terminate = b% accretor_overflow_terminate
         terminate_if_initial_overflow = b% terminate_if_initial_overflow
         terminate_if_L2_overflow = b% terminate_if_L2_overflow

         ! mass transfer controls
         mass_transfer_alpha = b% mass_transfer_alpha
         mass_transfer_beta = b% mass_transfer_beta
         mass_transfer_delta = b% mass_transfer_delta
         mass_transfer_gamma = b% mass_transfer_gamma
         limit_retention_by_mdot_edd = b% limit_retention_by_mdot_edd
         use_es_opacity_for_mdot_edd = b% use_es_opacity_for_mdot_edd
         use_this_for_mdot_edd_eta = b% use_this_for_mdot_edd_eta
         use_radiation_corrected_transfer_rate = b% use_radiation_corrected_transfer_rate
         initial_bh_spin = b% initial_bh_spin
         use_this_for_mdot_edd = b% use_this_for_mdot_edd
         mdot_scheme = b% mdot_scheme
         cur_mdot_frac = b% cur_mdot_frac
         max_explicit_abs_mdot = b% max_explicit_abs_mdot
         max_tries_to_achieve = b% max_tries_to_achieve
         solver_type = b% solver_type
         implicit_scheme_tolerance = b% implicit_scheme_tolerance
         implicit_scheme_tiny_factor = b% implicit_scheme_tiny_factor
         initial_change_factor = b% initial_change_factor
         change_factor_fraction = b% change_factor_fraction
         implicit_lambda = b% implicit_lambda
         max_change_factor = b% max_change_factor
         min_change_factor = b% min_change_factor
         num_tries_for_increase_change_factor = b% num_tries_for_increase_change_factor
         change_factor_increase = b% change_factor_increase
         starting_mdot = b% starting_mdot
         roche_min_mdot = b% roche_min_mdot
         min_mdot_for_implicit = b% min_mdot_for_implicit
         max_implicit_abs_mdot = b% max_implicit_abs_mdot
         report_rlo_solver_progress = b% report_rlo_solver_progress
         do_enhance_wind_1 = b% do_enhance_wind_1
         do_enhance_wind_2 = b% do_enhance_wind_2
         tout_B_wind_1 = b% tout_B_wind_1
         tout_B_wind_2 = b% tout_B_wind_2
         do_wind_mass_transfer_1 = b% do_wind_mass_transfer_1
         do_wind_mass_transfer_2 = b% do_wind_mass_transfer_2
         wind_BH_alpha_1 = b% wind_BH_alpha_1
         wind_BH_alpha_2 = b% wind_BH_alpha_2
         wind_BH_beta_1 = b% wind_BH_beta_1
         wind_BH_beta_2 = b% wind_BH_beta_2
         max_wind_transfer_fraction_1 = b% max_wind_transfer_fraction_1
         max_wind_transfer_fraction_2 = b% max_wind_transfer_fraction_2

         ! orbital jdot controls
         do_jdot_gr = b% do_jdot_gr
         do_jdot_ml = b% do_jdot_ml
         do_jdot_ls = b% do_jdot_ls
         do_jdot_missing_wind = b% do_jdot_missing_wind
         do_jdot_mb = b% do_jdot_mb
         include_accretor_mb = b% include_accretor_mb
         magnetic_braking_gamma = b% magnetic_braking_gamma
         keep_mb_on = b% keep_mb_on
         jdot_mb_min_qconv_env = b% jdot_mb_min_qconv_env
         jdot_mb_max_qconv_env = b% jdot_mb_max_qconv_env
         jdot_mb_max_qconv_core = b% jdot_mb_max_qconv_core
         jdot_mb_qlim_for_check_rad_core = b% jdot_mb_qlim_for_check_rad_core
         jdot_mb_qlim_for_check_conv_env = b% jdot_mb_qlim_for_check_conv_env
         jdot_mb_scale_for_low_qconv_env = b% jdot_mb_scale_for_low_qconv_env
         jdot_mb_mass_frac_for_scale = b% jdot_mb_mass_frac_for_scale
         jdot_multiplier = b% jdot_multiplier

         ! rotation and sync controls
         do_j_accretion = b% do_j_accretion
         do_tidal_sync = b% do_tidal_sync
         sync_type_1 = b% sync_type_1
         sync_type_2 = b% sync_type_2
         sync_mode_1 = b% sync_mode_1
         sync_mode_2 = b% sync_mode_2
         Ftid_1 = b% Ftid_1
         Ftid_2 = b% Ftid_2
         do_initial_orbit_sync_1 = b% do_initial_orbit_sync_1
         do_initial_orbit_sync_2 = b% do_initial_orbit_sync_2
         tidal_reduction = b% tidal_reduction
         
         ! eccentricity controls
         do_tidal_circ = b% do_tidal_circ
         circ_type_1 = b% circ_type_1
         circ_type_2 = b% circ_type_2
         use_eccentricity_enhancement = b% use_eccentricity_enhancement
         max_abs_edot_tidal = b% max_abs_edot_tidal
         max_abs_edot_enhance = b% max_abs_edot_enhance
         min_eccentricity = b% min_eccentricity
         max_eccentricity = b% max_eccentricity
         anomaly_steps = b% anomaly_steps

         ! irradiation controls
         accretion_powered_irradiation = b% accretion_powered_irradiation
         use_accretor_luminosity_for_irrad = b% use_accretor_luminosity_for_irrad
         col_depth_for_eps_extra = b% col_depth_for_eps_extra
         irrad_flux_at_std_distance = b% irrad_flux_at_std_distance
         std_distance_for_irradiation = b% std_distance_for_irradiation
         max_F_irr = b% max_F_irr

         !common envelope controls
         CE_alpha = b% CE_alpha
         CE_alpha_th = b% CE_alpha_th
         CE_alpha_core = b% CE_alpha_core
         CE_mass_loss_rate_high = b% CE_mass_loss_rate_high
         CE_mass_loss_rate_low = b% CE_mass_loss_rate_low
         CE_rel_rlo_for_detachment = b% CE_rel_rlo_for_detachment
         CE_years_detached_to_terminate = b% CE_years_detached_to_terminate
         CE_begin_at_max_implicit_abs_mdot = b% CE_begin_at_max_implicit_abs_mdot
         CE_xa_diff_to_terminate = b% CE_xa_diff_to_terminate
         CE_terminate_when_core_overflows = b% CE_terminate_when_core_overflows
         CE_min_period_in_minutes = b% CE_min_period_in_minutes
         CE_fixed_lambda = b% CE_fixed_lambda
         
         ! miscellaneous controls
         keep_donor_fixed = b% keep_donor_fixed
         mdot_limit_donor_switch = b% mdot_limit_donor_switch
         use_other_rlo_mdot = b% use_other_rlo_mdot
         use_other_check_implicit_rlo = b% use_other_check_implicit_rlo
         use_other_implicit_function_to_solve = b% use_other_implicit_function_to_solve
         use_other_tsync = b% use_other_tsync
         use_other_sync_spin_to_orbit = b% use_other_sync_spin_to_orbit
         use_other_mdot_edd = b% use_other_mdot_edd
         use_other_adjust_mdots = b% use_other_adjust_mdots
         use_other_accreted_material_j = b% use_other_accreted_material_j
         use_other_jdot_gr = b% use_other_jdot_gr
         use_other_jdot_ml = b% use_other_jdot_ml
         use_other_jdot_ls = b% use_other_jdot_ls
         use_other_jdot_missing_wind = b% use_other_jdot_missing_wind
         use_other_jdot_mb = b% use_other_jdot_mb
         use_other_extra_jdot = b% use_other_extra_jdot
         use_other_binary_wind_transfer = b% use_other_binary_wind_transfer
         use_other_edot_tidal = b% use_other_edot_tidal
         use_other_edot_enhance = b% use_other_edot_enhance
         use_other_extra_edot = b% use_other_extra_edot
         use_other_CE_init = b% use_other_CE_init
         use_other_CE_rlo_mdot = b% use_other_CE_rlo_mdot
         use_other_CE_binary_evolve_step = b% use_other_CE_binary_evolve_step
         use_other_CE_binary_finish_step = b% use_other_CE_binary_finish_step
         
      end subroutine set_binary_controls_for_writing
      
      subroutine write_binary_controls(io,ierr)
         integer, intent(in) :: io
         integer, intent(out) :: ierr
         write(io, nml=binary_controls, iostat=ierr)  
      end subroutine write_binary_controls


      subroutine get_binary_control(b, name, val, ierr)
         use utils_lib, only: StrUpCase
         type (binary_info), pointer :: b
         character(len=*),intent(in) :: name
         character(len=*), intent(out) :: val
         integer, intent(out) :: ierr
   
         character(len(name)) :: upper_name
         character(len=512) :: str
         integer :: iounit,iostat,ind,i
   
   
         ! First save current controls
         call set_binary_controls_for_writing(b, ierr)
         if(ierr/=0) return
   
         ! Write namelist to temporay file
         open(newunit=iounit,status='scratch')
         write(iounit,nml=binary_controls)
         rewind(iounit)
   
         ! Namelists get written in captials
         upper_name = StrUpCase(name)
         val = ''
         ! Search for name inside namelist
         do 
            read(iounit,'(A)',iostat=iostat) str
            ind = index(str,trim(upper_name))
            if( ind /= 0 ) then
               val = str(ind+len_trim(upper_name)+1:len_trim(str)-1) ! Remove final comma and starting =
               do i=1,len(val)
                  if(val(i:i)=='"') val(i:i) = ' '
               end do
               exit
            end if
            if(is_iostat_end(iostat)) exit
         end do   
   
         if(len_trim(val) == 0 .and. ind==0 ) ierr = -1
   
         close(iounit)
   
      end subroutine get_binary_control
   
      subroutine set_binary_control(b, name, val, ierr)
         type (binary_info), pointer :: b
         character(len=*), intent(in) :: name, val
         character(len=len(name)+len(val)+19) :: tmp
         integer, intent(out) :: ierr
   
         ! First save current controls
         call set_binary_controls_for_writing(b, ierr)
         if(ierr/=0) return
   
         tmp=''
         tmp = '&binary_controls '//trim(name)//'='//trim(val)//' /'
   
         ! Load into namelist
         read(tmp, nml=binary_controls)
   
         ! Add to star
         call store_binary_controls(b, ierr)
         if(ierr/=0) return
   
      end subroutine set_binary_control

      end module binary_ctrls_io

