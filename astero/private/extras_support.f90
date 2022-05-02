! ***********************************************************************
!
!   Copyright (C) 2010  The MESA Team
!
!   this file is part of mesa.
!
!   mesa is free software; you can redistribute it and/or modify
!   it under the terms of the gnu general library public license as published
!   by the free software foundation; either version 2 of the license, or
!   (at your option) any later version.
!
!   mesa is distributed in the hope that it will be useful, 
!   but without any warranty; without even the implied warranty of
!   merchantability or fitness for a particular purpose.  see the
!   gnu library general public license for more details.
!
!   you should have received a copy of the gnu library general public license
!   along with this software; if not, write to the free software
!   foundation, inc., 59 temple place, suite 330, boston, ma 02111-1307 usa
!
! ***********************************************************************
 
      module extras_support

      use star_lib
      use star_def
      use const_def
      use utils_lib
      use astero_support
      use astero_def
      use astero_def

      implicit none


      contains
      
      
      subroutine get_all_el_info(s, ierr)
         type (star_info), pointer :: s
         integer, intent(out) :: ierr
         logical :: store_model
         integer :: l
         store_model = .true.
         ierr = 0

         do l = 0, 3
            if (nl(l) > 0) then
               call get_one_el_info(s, l, &
                  nu_lower_factor*freq_target(l,1), &
                  nu_upper_factor*freq_target(l,nl(l)), &
                  iscan_factor(l)*nl(l), 1, nl(l), store_model, &
                  oscillation_code, ierr)
               if (ierr /= 0) return
               store_model = .false.
            end if
         end do

      end subroutine get_all_el_info
      

      integer function do_astero_extras_check_model(s, id)

         type (star_info), pointer :: s
         integer, intent(in) :: id
         
         integer :: max_el_for_chi2, ierr, i, j, l, n
         logical :: store_model, checking_age
         real(dp) :: age_limit, model_limit, err, target_l0, X, Y, Z, &
            frac, surface_X, surface_Z, chi2_freq_and_ratios_fraction, &
            remaining_years, prev_max_years, min_max
         
         include 'formats'
         
         do_astero_extras_check_model = keep_going
         astero_max_dt_next = 1d99
         chi2 = -1
         chi2_seismo = -1
         chi2_spectro = -1
         FeH = -1
         delta_nu_model = -1
         nu_max_model = -1
         a_div_r = -1
         correction_r = -1
         checking_age = &
            eval_chi2_at_target_age_only .or. include_age_in_chi2_spectro
         
         if (checking_age) then
            if (num_smaller_steps_before_age_target <= 0 .or. &
                dt_for_smaller_steps_before_age_target <= 0) then
               write(*,*) 'ERROR: must set num_smaller_steps_before_age_target'
               write(*,*) 'and dt_for_smaller_steps_before_age_target'
               call mesa_error(__FILE__,__LINE__)
            end if
            if (age_target > s% star_age) then
               remaining_years = age_target - s% star_age
               if (s% astero_using_revised_max_yr_dt) &
                  s% max_years_for_timestep = s% astero_revised_max_yr_dt
               n = floor(remaining_years/s% max_years_for_timestep + 1d-6)
               j = num_smaller_steps_before_age_target
               if (remaining_years <= s% max_years_for_timestep) then
                  write(*,'(a40,i6,f20.10)') '(age_target - star_age)/age_sigma', &
                     s% model_number, (age_target - s% star_age)/age_sigma
                  s% max_years_for_timestep = remaining_years
                  s% astero_using_revised_max_yr_dt = .true.
                  s% astero_revised_max_yr_dt = s% max_years_for_timestep
                  astero_max_dt_next = s% max_years_for_timestep*secyer
               else if (n <= j) then
                  write(*,'(a40,i6,f20.10)') '(age_target - star_age)/age_sigma', &
                     s% model_number, (age_target - s% star_age)/age_sigma
                  prev_max_years = s% max_years_for_timestep
                  i = floor(remaining_years/dt_for_smaller_steps_before_age_target + 1d-6)
                  if ((i+1d-9)*dt_for_smaller_steps_before_age_target < remaining_years) then
                     s% max_years_for_timestep = remaining_years/(i+1)
                  else
                     s% max_years_for_timestep = remaining_years/i
                  end if
                  min_max = prev_max_years*s% reduction_factor_for_max_timestep
                  if (s% max_years_for_timestep < min_max) &
                     s% max_years_for_timestep = min_max
                  if (.not. s% astero_using_revised_max_yr_dt) then
                     s% astero_using_revised_max_yr_dt = .true.
                     write(*,2) 'begin reducing max timestep prior to age target', &
                        s% model_number, remaining_years
                  else if (s% astero_revised_max_yr_dt > s% max_years_for_timestep) then
                     write(*,2) 'reducing max timestep prior to age target', &
                        s% model_number, remaining_years
                  else if (s% max_years_for_timestep <= dt_for_smaller_steps_before_age_target) then
                     i = floor(remaining_years/s% max_years_for_timestep + 1d-6)
                     write(*,3) 'remaining steps and years until age target', &
                        s% model_number, i, remaining_years
                  else 
                     write(*,2) 'remaining_years until age target', &
                        s% model_number, remaining_years
                  end if
                  s% astero_revised_max_yr_dt = s% max_years_for_timestep
                  if (s% dt_next/secyer > s% max_years_for_timestep) &
                     astero_max_dt_next = s% max_years_for_timestep*secyer
               end if
            else if (include_age_in_chi2_spectro) then
               write(*,'(a40,i6,f20.10)') '(age_target - star_age)/age_sigma', &
                  s% model_number, (age_target - s% star_age)/age_sigma
               if (abs(s% max_years_for_timestep - dt_for_smaller_steps_before_age_target) > &
                     dt_for_smaller_steps_before_age_target*1d-2) then
                  write(*,1) 'dt_for_smaller_steps_before_age_target', &
                     dt_for_smaller_steps_before_age_target
                  write(*,1) 'max_years_for_timestep', &
                     s% max_years_for_timestep
                  call mesa_error(__FILE__,__LINE__,'bad max_years_for_timestep')
               end if
            end if
         else
            if (s% star_age < min_age_for_chi2) return
            s% max_years_for_timestep = max_yrs_dt_when_cold
         end if

         if (include_age_in_chi2_spectro .and. s% star_age < min_age_for_chi2) return         
         if (eval_chi2_at_target_age_only .and. s% star_age < age_target) return
         
         delta_nu_model = s% delta_nu
         nu_max_model = s% nu_max
         
         chi2_seismo_delta_nu_fraction = &
            min(1d0, max(0d0, chi2_seismo_delta_nu_fraction))
         chi2_seismo_nu_max_fraction = &
            min(1d0, max(0d0, chi2_seismo_nu_max_fraction))
         chi2_seismo_r_010_fraction = &
            min(1d0, max(0d0, chi2_seismo_r_010_fraction))
         chi2_seismo_r_02_fraction = &
            min(1d0, max(0d0, chi2_seismo_r_02_fraction))
         chi2_seismo_freq_fraction = min(1d0, max(0d0, 1d0 - &
            (chi2_seismo_r_010_fraction + &
               chi2_seismo_r_02_fraction + &
               chi2_seismo_delta_nu_fraction + &
               chi2_seismo_nu_max_fraction)))
         
         if (s% L_nuc_burn_total < s% L_phot*Lnuc_div_L_limit .or. &
               s% star_age < min_age_limit) then
            return
         end if
         
         if (.not. checking_age) then
         
            age_limit = avg_age_top_samples + avg_age_sigma_limit*avg_age_sigma
            if (s% star_age > age_limit) then
               write(*,1) 'star age > limit from top samples', s% star_age, age_limit
               write(*,1) 'avg_age_top_samples', avg_age_top_samples
               write(*,1) 'avg_age_sigma_limit', avg_age_sigma_limit
               write(*,1) 'avg_age_sigma', avg_age_sigma
               do_astero_extras_check_model = terminate
               return
            end if
         
            model_limit = &
               avg_model_number_top_samples + &
                  avg_model_number_sigma_limit*avg_model_number_sigma
            if (dble(s% model_number) > model_limit) then
               write(*,2) 'model number > limit from top samples', &
                  s% model_number, model_limit
               write(*,1) 'avg_model_number_top_samples', avg_model_number_top_samples
               write(*,1) 'avg_model_number_sigma_limit', avg_model_number_sigma_limit
               write(*,1) 'avg_model_number_sigma', avg_model_number_sigma
               do_astero_extras_check_model = terminate
               return
            end if

         end if
         
         surface_X = max(s% surface_h1, 1d-10)
         surface_He = s% surface_he3 + s% surface_he4
         surface_Z = max(1d-99, min(1d0, 1 - (surface_X + surface_He)))
         surface_Z_div_X = surface_Z/surface_X
         FeH = log10((surface_Z_div_X)/Z_div_X_solar)
         logg = log10(s% grav(1))
         logR = log10(s% photosphere_r)
         if (.not. include_Rcz_in_chi2_spectro) then
            Rcz = 0
         else
            do i = 1, s% nz-1 ! locate bottom of solar convective zone
               if (s% mixing_type(i+1) /= convective_mixing &
                     .and. s% mixing_type(i) == convective_mixing) then
                  if (s% r(i+1) > 0.25*Rsun .and. s% r(i) < 0.9*Rsun) then
                     Rcz = s% r(i)/Rsun
                     exit
                  end if
               end if
            end do
         end if

         ! must set constraint values before checking limits
         do i = 1, max_constraints
            if (my_var_name(i) == '') cycle
            call star_astero_procs% set_my_vars(id, my_var_name(i), my_var(i), ierr)
            if (ierr /= 0) call mesa_error(__FILE__, __LINE__, 'ierr /=0 in set_my_vars')
         end do
         
         call check_limits
         if (do_astero_extras_check_model /= keep_going) return

         chi2_spectro = get_chi2_spectro(s)         
         if (is_bad(chi2_spectro)) then
            write(*,1) 'bad chi2_spectro', chi2_spectro
            write(*,1) 'FeH', FeH
            write(*,1) 'surface_Z', surface_Z
            write(*,1) 'surface_X', surface_X
            write(*,1) 'Z_div_X_solar', Z_div_X_solar
            chi2_spectro = 1d99
            do_astero_extras_check_model = terminate
            return
            !stop
         end if
         
         have_radial = .false.
         have_nonradial = .false.
         model_ratios_n = 0

         do l = 0, 3
            model_freq(l,1:nl(l)) = 0
         end do
         
         if (delta_nu_sigma > 0) then
            chi2_delta_nu = pow2((delta_nu - delta_nu_model)/delta_nu_sigma)
            if (trace_chi2_seismo_delta_nu_info) &
               write(*,1) 'chi2_delta_nu', chi2_delta_nu
         else
            chi2_delta_nu = 0
         end if
         
         chi2_nu_max = 0
         if (chi2_seismo_nu_max_fraction > 0) then
            if (nu_max <= 0) then
               write(*,2) 'must supply nu_max'
               do_astero_extras_check_model = terminate
               return
            end if         
            if (nu_max_sigma <= 0) then
               write(*,2) 'must supply nu_max_sigma'
               do_astero_extras_check_model = terminate
               return
            end if
            chi2_nu_max = pow2((nu_max - nu_max_model)/nu_max_sigma)
            if (trace_chi2_seismo_nu_max_info) &
               write(*,1) 'chi2_nu_max', chi2_nu_max 
         end if
         
         chi2_freq_and_ratios_fraction = &
            chi2_seismo_freq_fraction + &
            chi2_seismo_r_010_fraction + &
            chi2_seismo_r_02_fraction
            
         if (chi2_seismo_fraction <= 0d0) then
            ! no need to get frequencies
            chi2_seismo = &
               chi2_seismo_delta_nu_fraction*chi2_delta_nu + &
               chi2_seismo_nu_max_fraction*chi2_nu_max
            frac = chi2_seismo_fraction
            chi2 = frac*chi2_seismo + (1-frac)*chi2_spectro         
            write(*,'(a50,i6,99f16.2)') 'chi^2 combined, chi^2 seismo, chi^2 spectro', &
               s% model_number, chi2, chi2_seismo, chi2_spectro
            if (best_chi2 < 0 .or. chi2 < best_chi2) call save_best_info(s)
            if (chi2 < chi2_radial_limit .and. .not. checking_age) &
               s% max_years_for_timestep = max_yrs_dt_when_warm
            call final_checks
            return
         end if

         if (chi2_spectro > chi2_spectroscopic_limit) then
            write(*,'(a50,i6,99f16.2)') 'chi2_spectro > limit', &
                  s% model_number, chi2_spectro, chi2_spectroscopic_limit
            call check_too_many_bad
            return
         end if
         
         if (chi2_delta_nu > chi2_delta_nu_limit) then
            write(*,'(a50,i6,99f16.2)') 'chi2_delta_nu > limit', &
                  s% model_number, chi2_delta_nu, chi2_delta_nu_limit, &
                  delta_nu_model, delta_nu
            call check_too_many_bad
            return
         end if

         ! chi2_spectro <= limit and chi2_delta_nu <= limit        

         if (.not. checking_age) then
            s% max_years_for_timestep = max_yrs_dt_when_warm
            if (s% dt > max_yrs_dt_when_warm*secyer) then
               s% dt = max_yrs_dt_when_warm*secyer
               s% timestep_hold = s% model_number + 10
               write(*,'(a50,i6,1p,99e16.4)') &
                  'redo with smaller timestep for "warm" limit', &
                     s% model_number, max_yrs_dt_when_warm
               do_astero_extras_check_model = redo
               return
            end if            
         end if
         
         if (.not. checking_age) then
            s% max_years_for_timestep = max_yrs_dt_when_hot 
            if (s% dt > max_yrs_dt_when_hot*secyer) then
               s% dt = max_yrs_dt_when_hot*secyer
               s% timestep_hold = s% model_number + 10
               write(*,'(a50,i6,1p,99e16.4)') &
                  'redo with smaller timestep for "hot" limit', &
                  s% model_number, max_yrs_dt_when_hot
               do_astero_extras_check_model = redo
               return
            end if
         end if
         
         store_model = .true.
         if (nl(0) > 0 .and. chi2_freq_and_ratios_fraction > 0d0) then
            if (.not. get_radial(oscillation_code)) then
               !write(*,'(a65,i6)') 'failed to find all required l=0 modes', s% model_number
               if (trace_chi2_seismo_frequencies_info) then
                  write(*,2) 'results for l=0'
                  i = 0
                  do j = 1, num_results
                     if (el(j) /= 0) cycle
                     i = i+1
                     write(*,2) 'freq', i, cyclic_freq(j)
                  end do
                  write(*,'(A)')
               end if
               call check_too_many_bad
               return 
            end if
            store_model = .false.
            have_radial = .true.
         end if

         do l = 1, 3
            !write(*,3) 'l, nl(l)', l, nl(l)
            if (nl(l) > 0 .and. chi2_freq_and_ratios_fraction > 0d0) then
               call get_one_el_info(s, l, &
                  nu_lower_factor*freq_target(l,1), &
                  nu_upper_factor*freq_target(l,nl(l)), &
                  iscan_factor(l)*nl(l), 1, nl(l), store_model, &
                  oscillation_code, ierr)
               if (ierr /= 0) then
                  if (trace_chi2_seismo_frequencies_info) write(*,*)
                  write(*,'(a65,i4,i6)') 'failed to find all required modes', l, s% model_number
                  if (trace_chi2_seismo_frequencies_info) then
                     write(*,2) 'results for l =', l
                     i = 0
                     do j = 1, num_results
                        if (el(j) /= l) cycle
                        i = i+1
                        write(*,2) 'freq', i, cyclic_freq(j)
                     end do
                     write(*,'(A)')
                  end if
                  call check_too_many_bad
                  return
               end if
               store_model = .false.
            end if
         end do
         
         have_nonradial = .true.
         
         if (chi2_freq_and_ratios_fraction > 0d0) then
            call get_freq_corr(s, .false., ierr)
            if (ierr /= 0) then
               write(*,'(a65,i6)') 'failed in get_freq_corr', s% model_number
               return
            end if
            if (nl(3) > 0) then
               max_el_for_chi2 = 3
            else if (nl(2) > 0) then
               max_el_for_chi2 = 2
            else if (nl(1) > 0) then
               max_el_for_chi2 = 1
            else
               max_el_for_chi2 = 0
            end if
         else
            max_el_for_chi2 = -1
         end if
         
         !write(*,2) 'max_el_for_chi2', max_el_for_chi2
         if (chi2_seismo_r_010_fraction > 0 .and. max_el_for_chi2 >= 1) then

            call get_frequency_ratios( &
               .false., nl(0), model_freq_corr(0,:), nl(1), model_freq_corr(1,:), &
               model_ratios_n, model_ratios_l0_first, model_ratios_l1_first, &
               model_ratios_r01, model_ratios_r10)
            
            if (model_ratios_n /= ratios_n .or. &
                  model_ratios_l0_first /= ratios_l0_first .or. &
                     model_ratios_l1_first /= ratios_l1_first) then
               ierr = -1
               write(*,'(a,i6)') 'cannot calculate necessary chi^2 ratios for this model', s% model_number
               if (model_ratios_n /= ratios_n) &
                  write(*,*) '              model_ratios_n /= ratios_n', model_ratios_n, ratios_n
               if (model_ratios_l0_first /= ratios_l0_first) &
                  write(*,*) 'model_ratios_l0_first /= ratios_l0_first', &
                     model_ratios_l0_first, ratios_l0_first
               if (model_ratios_l1_first /= ratios_l1_first) &
                  write(*,*) 'model_ratios_l1_first /= ratios_l1_first', &
                     model_ratios_l1_first, ratios_l1_first
               call check_too_many_bad
               return               
            end if
            
         end if
               
         if (chi2_seismo_r_02_fraction > 0 .and. max_el_for_chi2 >= 2) then
            call get_r02_frequency_ratios( &
               .false., nl(0), model_freq_corr(0,:), nl(1), model_freq_corr(1,:), nl(2), model_freq_corr(2,:), model_ratios_r02)
         end if
         
         chi2 = get_chi2(s, max_el_for_chi2, .true., ierr)
         if (ierr /= 0) then
            write(*,'(a40,i6)') 'failed to calculate chi^2', s% model_number
            call check_too_many_bad
            return
         end if
         write(*,'(a50,i6,99f16.2)') 'chi^2 total, chi^2 radial', &
            s% model_number, chi2, chi2_radial
         
         if (use_other_after_get_chi2) then
            ierr = 0
            call astero_other_procs% other_after_get_chi2(s% id, ierr)
            if (ierr /= 0) then
               do_astero_extras_check_model = terminate
               return
            end if
         end if
          
         if (checking_age) then
            ! leave max_years_for_timestep as is
         else if (chi2 <= chi2_limit_for_smallest_timesteps) then
            s% max_years_for_timestep = max_yrs_dt_chi2_smallest_limit 
            if (s% dt > max_yrs_dt_chi2_smallest_limit*secyer) then
               s% dt = max_yrs_dt_chi2_smallest_limit*secyer
               s% timestep_hold = s% model_number + 10
               write(*,'(a50,i6,2f16.2,1p,99e16.4)') 'redo timestep for "smallest" chi^2 limit', &
                  s% model_number, chi2, chi2_limit_for_smallest_timesteps, &
                  max_yrs_dt_chi2_smallest_limit
               do_astero_extras_check_model = redo
               return
            end if         
         else if (chi2 <= chi2_limit_for_smaller_timesteps) then
            s% max_years_for_timestep = max_yrs_dt_chi2_smaller_limit 
            if (s% dt > max_yrs_dt_chi2_smaller_limit*secyer) then
               s% dt = max_yrs_dt_chi2_smaller_limit*secyer
               s% timestep_hold = s% model_number + 10
               write(*,'(a50,i6,2f16.2,1p,99e16.4)') 'redo timestep for "smaller" chi^2 limit', &
                  s% model_number, chi2, chi2_limit_for_smaller_timesteps, &
                  max_yrs_dt_chi2_smaller_limit
               do_astero_extras_check_model = redo
               return
            end if         
         else if (chi2 <= chi2_limit_for_small_timesteps) then
            s% max_years_for_timestep = max_yrs_dt_chi2_small_limit 
            if (s% dt > max_yrs_dt_chi2_small_limit*secyer) then
               s% dt = max_yrs_dt_chi2_small_limit*secyer
               s% timestep_hold = s% model_number + 10
               write(*,'(a50,i6,2f16.2,1p,99e16.4)') 'redo timestep for "small" chi^2 limit', &
                  s% model_number, chi2, chi2_limit_for_small_timesteps, &
                  max_yrs_dt_chi2_small_limit
               do_astero_extras_check_model = redo
               return
            end if         
         end if
         
         if (best_chi2 <= 0 .or. chi2 < best_chi2) then
            call save_best_info(s)
         end if
            
         call final_checks

         
         contains
         
         
         subroutine check_too_many_bad
            if (best_chi2 > 0) then
               num_chi2_too_big = num_chi2_too_big + 1
               if (num_chi2_too_big > limit_num_chi2_too_big) then
                  write(*,*) 'have reached too many bad chi2 limit'
                  do_astero_extras_check_model = terminate
               end if
               return
            end if
            num_chi2_too_big = 0
         end subroutine check_too_many_bad
         
         
         subroutine final_checks
            if (include_age_in_chi2_spectro .and. s% star_age >= max_age_for_chi2) then
               write(*,*) 'have reached max_age_for_chi2'
               do_astero_extras_check_model = terminate
            end if
            if (eval_chi2_at_target_age_only .and. s% star_age >= age_target) then
               write(*,*) 'have reached age_target'
               do_astero_extras_check_model = terminate
            end if
            if (best_chi2 > 0) then
               if (best_chi2 <= chi2_search_limit1 .and. &
                  chi2 >= chi2_search_limit2) then
                  write(*,*) 'have reached chi2_search_limit2'
                  do_astero_extras_check_model = terminate
                  return
               end if
               if (chi2 >= chi2_relative_increase_limit*best_chi2) then
                  num_chi2_too_big = num_chi2_too_big + 1
                  if (num_chi2_too_big > limit_num_chi2_too_big) then
                     write(*,*) 'have reached too many bad chi2 limit'
                     do_astero_extras_check_model = terminate
                  end if
                  return
               end if
               num_chi2_too_big = 0
            end if
         end subroutine final_checks
         

         logical function get_radial(code)
            character (len=*), intent(in) :: code
            integer :: ierr
            real(dp) :: chi2term
            include 'formats'
            ierr = 0
            get_radial = .false.
            model_freq(0,:) = 0d0
            call get_one_el_info(s, 0, &
               nu_lower_factor*freq_target(0,1), &
               nu_upper_factor*freq_target(0,nl(0)), &
               iscan_factor(0)*nl(0), 1, nl(0), store_model, &
               code, ierr)
            if (ierr /= 0) then
               !write(*,'(a65,i6)') 'failed in oscillation code', s% model_number
               !stop
               return
            end if
            if (.not. have_all_l0_freqs()) then
               write(*,'(a65,i6)') 'failed to find all required l=0 frequencies', &
                  s% model_number
               model_freq_corr(0,:) = model_freq(0,:) ! for plotting
               return
            end if
            call get_freq_corr(s, .true., ierr)
            if (ierr /= 0) then
               write(*,'(a65,i6)') 'failed in get_freq_corr', s% model_number
               return
            end if
            ! you might think we can use the function get_chi2 from astero_support for this
            ! but that causes a lot of problems
            ! chi2_radial = get_chi2(s, 0, .false., ierr)
            ! instead, we just add up the radial chi^2
            chi2_radial = 0.0_dp
            n = 0
            do i = 1, nl(0)
               if (freq_target(0,i) < 0) cycle
               chi2term = &
                  pow2((model_freq_corr(0,i) - freq_target(0,i))/freq_sigma(0,i))
               if (.false. .and. trace_chi2_seismo_frequencies_info) &
                  write(*,'(4i6,99(1pe20.10))') &
                     s% model_number, i, 0, model_order(0,i), chi2term, model_freq(0,i), &
                     model_freq_corr(0,i), freq_target(0,i), freq_sigma(0,i), safe_log10(model_inertia(0,i))
               chi2_radial = chi2_radial + chi2term
               n = n + 1
            end do

            if (normalize_chi2_seismo_frequencies) then
               chi2_radial = chi2_radial/max(1,n)
            end if

            if (ierr /= 0) then
               write(*,'(a65,i6)') 'failed to get chi2_radial', s% model_number
               return
            end if
            !write(*,1) trim(code) // ' chi2_radial', chi2_radial
            if (chi2_radial > chi2_radial_limit .and. nl(1) + nl(2) + nl(3) > 0) then
               write(*,'(a50,i6,99f16.2)') &
                  'chi2_radial > chi2_radial_limit', &
                  s% model_number, chi2_radial, chi2_radial_limit
               return
            end if
            get_radial = .true.
         end function get_radial
         
         
         logical function have_all_l0_freqs()
            integer :: i, cnt
            real(dp) :: prev
            cnt = 0
            have_all_l0_freqs = .true.
            if (nl(0) <= 0) return
            prev = model_freq(0,1)
            do i=2,nl(0)
               if (freq_target(0,i) < 0) cycle
               if (model_freq(0,i) == prev) then
                  have_all_l0_freqs = .false.
                  if (cnt == 0) write(*,'(i30,4x,a)',advance='no') &
                     s% model_number, 'missing l=0 freq number:'
                  cnt = cnt+1
                  write(*,'(i3)',advance='no') i
               end if
            end do
            if (cnt > 0) write(*,*)
         end function have_all_l0_freqs
         
         
         subroutine check_limits
            real(dp) :: logg_limit, logL_limit, Teff_limit, delta_nu_limit, &
               logR_limit, surface_Z_div_X_limit, surface_He_limit, Rcz_limit, &
               my_var_limit
            integer :: nz, i
            include 'formats'
            nz = s% nz

            if (sigmas_coeff_for_Teff_limit /= 0 .and. Teff_sigma > 0) then
               Teff_limit = Teff_target + Teff_sigma*sigmas_coeff_for_Teff_limit
               if ((sigmas_coeff_for_Teff_limit > 0 .and. s% Teff > Teff_limit) .or. &
                   (sigmas_coeff_for_Teff_limit < 0 .and. s% Teff < Teff_limit)) then
                  write(*,*) 'have reached Teff limit'
                  write(*,1) 'Teff', s% Teff
                  write(*,1) 'Teff_limit', Teff_limit
                  write(*,'(A)')
                  do_astero_extras_check_model = terminate
                  return
               end if    
               if (trace_limits) then
                  write(*,1) 'Teff', s% Teff
                  write(*,1) 'Teff_limit', Teff_limit
               end if
            end if
            
            if (sigmas_coeff_for_logg_limit /= 0 .and. logg_sigma > 0) then    
               logg_limit = logg_target + logg_sigma*sigmas_coeff_for_logg_limit
               if ((sigmas_coeff_for_logg_limit > 0 .and. logg > logg_limit) .or. &
                   (sigmas_coeff_for_logg_limit < 0 .and. logg < logg_limit)) then
                  write(*,*) 'have reached logg limit'
                  write(*,1) 'logg', logg
                  write(*,1) 'logg_limit', logg_limit
                  write(*,'(A)')
                  do_astero_extras_check_model = terminate
                  return
               end if
               if (trace_limits) then
                  write(*,1) 'logg', logg
                  write(*,1) 'logg_limit', logg_limit
               end if
            end if
            
            if (sigmas_coeff_for_logL_limit /= 0 .and. logL_sigma > 0) then
               logL_limit = logL_target + logL_sigma*sigmas_coeff_for_logL_limit
               if ((sigmas_coeff_for_logL_limit > 0 .and. s% log_surface_luminosity > logL_limit) .or. &
                   (sigmas_coeff_for_logL_limit < 0 .and. s% log_surface_luminosity < logL_limit)) then
                  write(*,*) 'have reached logL limit'
                  write(*,1) 'logL', s% log_surface_luminosity
                  write(*,1) 'logL_limit', logL_limit
                  write(*,'(A)')
                  do_astero_extras_check_model = terminate
                  return
               end if
               if (trace_limits) then
                  write(*,1) 'logL', s% log_surface_luminosity
                  write(*,1) 'logL_limit', logL_limit
               end if
            end if
            
            if (sigmas_coeff_for_delta_nu_limit /= 0 .and. delta_nu_sigma > 0 .and. delta_nu > 0) then
               delta_nu_limit = &
                  delta_nu + delta_nu_sigma*sigmas_coeff_for_delta_nu_limit
               if ((sigmas_coeff_for_delta_nu_limit > 0 .and. delta_nu_model > delta_nu_limit) .or. &
                   (sigmas_coeff_for_delta_nu_limit < 0 .and. delta_nu_model < delta_nu_limit)) then
                  write(*,*) 'have reached delta_nu limit'
                  write(*,1) 'delta_nu_model', delta_nu_model
                  write(*,1) 'delta_nu_limit', delta_nu_limit
                  write(*,'(A)')
                  do_astero_extras_check_model = terminate
                  return
               end if
               if (trace_limits) then
                  write(*,1) 'delta_nu_model', delta_nu_model
                  write(*,1) 'delta_nu_limit', delta_nu_limit
               end if
            end if
            
            if (sigmas_coeff_for_logR_limit /= 0 .and. logR_sigma > 0) then
               logR_limit = logR_target + logR_sigma*sigmas_coeff_for_logR_limit
               if ((sigmas_coeff_for_logR_limit > 0 .and. logR > logR_limit) .or. &
                   (sigmas_coeff_for_logR_limit < 0 .and. logR < logR_limit)) then
                  write(*,*) 'have reached logR limit'
                  write(*,1) 'logR', logR
                  write(*,1) 'logR_limit', logR_limit
                  write(*,'(A)')
                  do_astero_extras_check_model = terminate
                  return
               end if
               if (trace_limits) then
                  write(*,1) 'logR', logR
                  write(*,1) 'logR_limit', logR_limit
               end if
            end if
            
            if (sigmas_coeff_for_surface_Z_div_X_limit /= 0 .and. surface_Z_div_X_sigma > 0) then
               surface_Z_div_X_limit = surface_Z_div_X_target + &
                  surface_Z_div_X_sigma*sigmas_coeff_for_surface_Z_div_X_limit
               if ((sigmas_coeff_for_surface_Z_div_X_limit > 0 .and. &
                     surface_Z_div_X > surface_Z_div_X_limit) .or. &
                     (sigmas_coeff_for_surface_Z_div_X_limit < 0 .and. &
                      surface_Z_div_X < surface_Z_div_X_limit)) then
                  write(*,*) 'have reached surface_Z_div_X limit'
                  write(*,1) 'surface_Z_div_X', surface_Z_div_X
                  write(*,1) 'surface_Z_div_X_limit', surface_Z_div_X_limit
                  write(*,'(A)')
                  do_astero_extras_check_model = terminate
                  return
               end if
               if (trace_limits) then
                  write(*,1) 'surface_Z_div_X', surface_Z_div_X
                  write(*,1) 'surface_Z_div_X_limit', surface_Z_div_X_limit
               end if
            end if
            
            if (sigmas_coeff_for_surface_He_limit /= 0 .and. surface_He_sigma > 0) then
               surface_He_limit = surface_He_target + &
                     surface_He_sigma*sigmas_coeff_for_surface_He_limit
               if ((sigmas_coeff_for_surface_He_limit > 0 .and. &
                     surface_He > surface_He_limit) .or. &
                   (sigmas_coeff_for_surface_He_limit < 0 .and. &
                     surface_He < surface_He_limit)) then
                  write(*,*) 'have reached surface_He limit'
                  write(*,1) 'surface_He', surface_He
                  write(*,1) 'surface_He_limit', surface_He_limit
                  write(*,'(A)')
                  do_astero_extras_check_model = terminate
                  return
               end if
               if (trace_limits) then
                  write(*,1) 'surface_He', surface_He
                  write(*,1) 'surface_He_limit', surface_He_limit
               end if
            end if
            
            if (sigmas_coeff_for_Rcz_limit /= 0 .and. Rcz_sigma > 0) then
               Rcz_limit = Rcz_target + Rcz_sigma*sigmas_coeff_for_Rcz_limit
               if ((sigmas_coeff_for_Rcz_limit > 0 .and. Rcz > Rcz_limit) .or. &
                   (sigmas_coeff_for_Rcz_limit < 0 .and. Rcz < Rcz_limit)) then
                  write(*,*) 'have reached Rcz limit'
                  write(*,1) 'Rcz', Rcz
                  write(*,1) 'Rcz_limit', Rcz_limit
                  write(*,'(A)')
                  do_astero_extras_check_model = terminate
                  return
               end if
               if (trace_limits) then
                  write(*,1) 'Rcz', Rcz
                  write(*,1) 'Rcz_limit', Rcz_limit
               end if
            end if

            do i = 1, max_constraints
               if (sigmas_coeff_for_my_var_limit(i) /= 0 .and. my_var_sigma(i) > 0) then
                  my_var_limit = &
                     my_var_target(i) + my_var_sigma(i)*sigmas_coeff_for_my_var_limit(i)
                  if ((sigmas_coeff_for_my_var_limit(i) > 0 .and. &
                           my_var(i) > my_var_limit) .or. &
                      (sigmas_coeff_for_my_var_limit(i) < 0 .and. &
                           my_var(i) < my_var_limit)) then
                     write(*,*) 'have reached my_var(i) limit'
                     write(*,1) 'my_var(i)', my_var(i)
                     write(*,1) 'my_var(i)_limit', my_var_limit
                     write(*,'(A)')
                     do_astero_extras_check_model = terminate
                     return
                  end if
                  if (trace_limits) then
                     write(*,1) 'my_var(i)', my_var(i)
                     write(*,1) 'my_var_limit', my_var_limit
                  end if
               end if
            end do
            
         end subroutine check_limits         

      end function do_astero_extras_check_model
      
      
      real(dp) function get_chi2_spectro(s)
         type (star_info), pointer :: s
         integer :: cnt, i
         real(dp) :: logL, sum
         include 'formats'
         cnt = 0
         sum = 0
         if (include_logL_in_chi2_spectro) then
            cnt = cnt + 1
            logL = s% log_surface_luminosity
            sum = sum + pow2((logL - logL_target)/logL_sigma)
         end if
         if (include_logg_in_chi2_spectro) then
            cnt = cnt + 1
            sum = sum + pow2((logg - logg_target)/logg_sigma)
         end if
         if (include_Teff_in_chi2_spectro) then
            cnt = cnt + 1
            sum = sum + pow2((s% Teff - Teff_target)/Teff_sigma)
         end if
         if (include_FeH_in_chi2_spectro) then
            cnt = cnt + 1
            sum = sum + pow2((FeH - FeH_target)/FeH_sigma)
         end if
         if (include_logR_in_chi2_spectro) then
            cnt = cnt + 1
            sum = sum + pow2((logR - logR_target)/logR_sigma)
         end if
         if (include_age_in_chi2_spectro) then
            cnt = cnt + 1
            sum = sum + pow2((s% star_age - age_target)/age_sigma)
         end if
         if (include_surface_Z_div_X_in_chi2_spectro) then
            cnt = cnt + 1
            sum = sum + &
               pow2( &
               (surface_Z_div_X - surface_Z_div_X_target)/surface_Z_div_X_sigma)
         end if
         if (include_surface_He_in_chi2_spectro) then
            cnt = cnt + 1
            sum = sum + pow2( &
               (surface_He - surface_He_target)/surface_He_sigma)
         end if
         if (include_Rcz_in_chi2_spectro) then
            cnt = cnt + 1
            sum = sum + pow2((Rcz - Rcz_target)/Rcz_sigma)
         end if

         do i = 1, max_constraints
            if (include_my_var_in_chi2_spectro(i)) then
               cnt = cnt + 1
               sum = sum + pow2( &
                  (my_var(i) - my_var_target(i))/my_var_sigma(i))
            end if
         end do

         if (normalize_chi2_spectro) then
            get_chi2_spectro = sum/cnt
         else
            get_chi2_spectro = sum
         end if
      end function get_chi2_spectro
      
      
      subroutine store_best_info(s)
         type (star_info), pointer :: s
         integer :: i, l
      
         best_chi2 = chi2
         best_chi2_seismo = chi2_seismo
         best_chi2_spectro = chi2_spectro
         
         best_age = s% star_age
         best_model_number = s% model_number
         best_radius = s% photosphere_r
         best_logL = s% log_surface_luminosity
         best_Teff = s% Teff
         best_logg = logg
         best_FeH = FeH
         
         best_logR = logR
         best_surface_Z_div_X = surface_Z_div_X
         best_surface_He = surface_He
         best_Rcz = Rcz
         best_my_var(1:max_constraints) = my_var(1:max_constraints)

         best_my_param1 = my_param1
         best_my_param2 = my_param2
         best_my_param3 = my_param3
         
         best_delta_nu = delta_nu_model
         best_nu_max = nu_max_model
         best_surf_coef1 = surf_coef1
         best_surf_coef2 = surf_coef2

         do l = 0, 3
            do i = 1, nl(l)
               best_order(l,i) = model_order(l,i)
               best_freq(l,i) = model_freq(l,i)
               best_freq_corr(l,i) = model_freq_corr(l,i)
               best_inertia(l,i) = model_inertia(l,i)
            end do
         end do
         
         best_ratios_r01(:) = 0d0
         best_ratios_r10(:) = 0d0
         best_ratios_r02(:) = 0d0

         do i=1,ratios_n
            best_ratios_r01(i) = model_ratios_r01(i)
            best_ratios_r10(i) = model_ratios_r10(i)
         end do
         
         do i=1,nl(0)
            best_ratios_r02(i) = model_ratios_r02(i)
         end do
      
      end subroutine store_best_info


      subroutine set_current_from_best(s)
         type (star_info), pointer :: s
         integer :: i, l
      
         chi2 = best_chi2
         chi2_seismo = best_chi2_seismo
         chi2_spectro = best_chi2_spectro
         
         delta_nu_model = best_delta_nu
         nu_max_model = best_nu_max
         surf_coef1 = best_surf_coef1
         surf_coef2 = best_surf_coef2

         do l = 0, 3
            do i = 1, nl(l)
               model_order(l,i) = best_order(l,i)
               model_freq(l,i) = best_freq(l,i)
               model_freq_corr(l,i) = best_freq_corr(l,i)
               model_inertia(l,i) = best_inertia(l,i)
            end do
         end do
         
         do i=1,ratios_n
            model_ratios_r01(i) = best_ratios_r01(i)
            model_ratios_r10(i) = best_ratios_r10(i)
         end do
         
         do i=1,nl(0)
            model_ratios_r02(i) = best_ratios_r02(i)
         end do
      
      end subroutine set_current_from_best
      
      
      subroutine save_best_info(s)     
         use pgstar_astero_plots, only: write_plot_to_file
         type (star_info), pointer :: s
         integer :: ierr
         logical :: write_controls_info_with_profile
         character (len=256) :: filename
         
         include 'formats'
         
         if (save_model_for_best_model) then
            ierr = 0
            filename = trim(astero_results_directory) // '/' // trim(best_model_save_model_filename)
            if (.not. folder_exists(trim(astero_results_directory))) call mkdir(trim(astero_results_directory))
            call star_write_model(s% id, filename, ierr)
            if (ierr /= 0) then
               write(*,*) 'failed in star_write_model'
               call mesa_error(__FILE__,__LINE__)
            end if
            write(*, '(a,i7)') 'save ' // filename, s% model_number
         end if
         
         if (write_fgong_for_best_model) then
            ierr = 0
            filename = trim(astero_results_directory) // '/' // trim(best_model_fgong_filename)
            if (.not. folder_exists(trim(astero_results_directory))) call mkdir(trim(astero_results_directory))
            call star_export_pulse_data(s%id, 'FGONG', filename, &
               add_center_point, keep_surface_point, add_atmosphere, ierr)
            if (ierr /= 0) then
               write(*,*) 'failed in star_export_pulse_data'
               call mesa_error(__FILE__,__LINE__)
            end if
         end if
         
         if (write_gyre_for_best_model) then
            ierr = 0
            filename = trim(astero_results_directory) // '/' // trim(best_model_gyre_filename)
            if (.not. folder_exists(trim(astero_results_directory))) call mkdir(trim(astero_results_directory))
            call star_export_pulse_data(s%id, 'GYRE', filename, &
               add_center_point, keep_surface_point, add_atmosphere, ierr)
            if (ierr /= 0) then
               write(*,*) 'failed in star_export_pulse_data'
               call mesa_error(__FILE__,__LINE__)
            end if
         end if
         
         if (write_profile_for_best_model) then
            ierr = 0
            filename = trim(astero_results_directory) // '/' // trim(best_model_profile_filename)
            if (.not. folder_exists(trim(astero_results_directory))) call mkdir(trim(astero_results_directory))
            write_controls_info_with_profile = s% write_controls_info_with_profile
            s% write_controls_info_with_profile = .false.
            call star_write_profile_info(s% id, filename, ierr)
            s% write_controls_info_with_profile = write_controls_info_with_profile
            if (ierr /= 0) then
               write(*,*) 'failed in star_write_profile_info'
               call mesa_error(__FILE__,__LINE__)
            end if
            call save_profile(s% id, 3, ierr)
            if (ierr /= 0) then
               write(*,*) 'failed in save_profile'
               call mesa_error(__FILE__,__LINE__)
            end if
         end if         
         
         if (len_trim(echelle_best_model_file_prefix) > 0) then
            ! note: sample_number hasn't been incremented yet so must add 1
            call write_plot_to_file( &
               s, p_echelle, echelle_best_model_file_prefix, sample_number+1, ierr)        
         end if
         
         if (len_trim(ratios_best_model_file_prefix) > 0) then
            ! note: sample_number hasn't been incremented yet so must add 1
            call write_plot_to_file( &
               s, p_ratios, ratios_best_model_file_prefix, sample_number+1, ierr)   
         end if
         
         call store_best_info(s)
            
      end subroutine save_best_info
      
      
      subroutine write_best(num)
         integer, intent(in) :: num
         integer :: ierr, iounit
         character (len=256) :: format_string, num_string, filename
         integer, parameter :: max_len_out = 2000
         character (len=max_len_out) :: script         
         ierr = 0
         iounit = alloc_iounit(ierr)
         if (ierr /= 0) return

         if (.not. folder_exists(trim(astero_results_directory))) call mkdir(trim(astero_results_directory))

         write(format_string,'( "(i",i2.2,".",i2.2,")" )') num_digits, num_digits
         write(num_string,format_string) num
         filename = trim(astero_results_directory) // '/' // trim(sample_results_prefix) // trim(num_string) // trim(sample_results_postfix)
         open(unit=iounit, file=trim(filename), action='write', status='replace', iostat=ierr)
         if (ierr == 0) then
            call show_best(iounit)
            close(iounit)
            write(*,*) 'save best model results to ' // trim(filename)
         else
            write(*,*) 'failed to open ' // trim(filename)
         end if
         call free_iounit(iounit)
      end subroutine write_best


      integer function astero_extras_check_model(id)
            
         integer, intent(in) :: id
         integer :: other_check, ierr
         type (star_info), pointer :: s
         
         include 'formats'
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         if (s% job% astero_just_call_my_extras_check_model) then
            astero_extras_check_model = star_astero_procs% extras_check_model(id)
            best_chi2 = 0
         else
            other_check = star_astero_procs% extras_check_model(id)
            astero_extras_check_model = &
                  do_astero_extras_check_model(s, id)
            if (other_check > astero_extras_check_model) &
               astero_extras_check_model = other_check
         end if
         
         star_model_number = s% model_number
         if (star_model_number /= save_mode_model_number) return
         call get_all_el_info(s,ierr)
               
      end function astero_extras_check_model

      
      integer function astero_extras_finish_step(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         if (s% job% pgstar_flag) then
            ierr = 0
            call read_astero_pgstar_controls(inlist_astero_fname, ierr)
            if (ierr /= 0) then
               astero_extras_finish_step = terminate
               return
            end if
         end if
         astero_extras_finish_step = star_astero_procs% extras_finish_step(id)
         call store_extra_info(s)
         
         s% dt_next = min(s% dt_next, astero_max_dt_next)
         
      end function astero_extras_finish_step
      

      subroutine astero_extras_controls(id, ierr)
         !use run_star_extras, only: extras_controls, will_set_my_param
         use pgstar_astero_plots, only: astero_pgstar_plots_info
         use gyre_support, only: gyre_is_enabled, init_gyre
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         real(dp) :: X, Y, Z, FeH, f_ov, a, b, c
         type (star_info), pointer :: s
         include 'formats'
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         
         write(*,*) 'enter astero_extras_controls'
         
         call star_astero_procs% extras_controls(id, ierr)
         if (ierr /= 0) return
         
         
         s% extras_startup => astero_extras_startup
         s% extras_check_model => astero_extras_check_model
         s% extras_finish_step => astero_extras_finish_step
         s% extras_after_evolve => astero_extras_after_evolve
         s% how_many_extra_history_columns => astero_how_many_extra_history_columns
         s% data_for_extra_history_columns => astero_data_for_extra_history_columns
         s% how_many_extra_profile_columns => astero_how_many_extra_profile_columns
         s% data_for_extra_profile_columns => astero_data_for_extra_profile_columns  
         
         
         if (s% job% astero_just_call_my_extras_check_model) return
         
         s% other_pgstar_plots_info => astero_pgstar_plots_info
         s% use_other_pgstar_plots = .true.
         
         
         ! overwrite various inlist controls

         if (vary_alpha) then
            s% mixing_length_alpha = next_alpha_to_try
         else
            s% mixing_length_alpha = first_alpha
         end if

         if (vary_f_ov) then
            f_ov = next_f_ov_to_try
         else
            f_ov = first_f_ov
         end if
      
         if (vary_FeH) then
            FeH = next_FeH_to_try
         else
            FeH = first_FeH
         end if

         initial_FeH = FeH
         initial_Z_div_X = Z_div_X_solar*exp10(FeH)

         if (Y_depends_on_Z) then
            a = initial_Z_div_X
            b = dYdZ
            c = 1d0 + a*(1d0 + b)
            X = (1d0 - Y0)/c
            Y = (Y0 + a*(b + Y0))/c
            Z = 1d0 - (X + Y)
            !write(*,1) 'init X', X
            !write(*,1) 'init Y', Y
            !write(*,1) 'init Z', Z
            !stop
         else 
            if (vary_Y) then
               Y = next_Y_to_try
            else
               Y = first_Y
            end if
            X = (1d0 - Y)/(1d0 + initial_Z_div_X)
            Z = X*initial_Z_div_X
         end if

         if (vary_mass) then
            s% job% new_mass = next_mass_to_try
         else
            s% job% new_mass = first_mass
         end if

         if (vary_my_param1) then
            call star_astero_procs% will_set_my_param( &
               s% id, 1, next_my_param1_to_try, ierr)
            if (ierr /= 0) return
            my_param1 = next_my_param1_to_try
         else
            call star_astero_procs% will_set_my_param( &
               s% id, 1, first_my_param1, ierr)
            if (ierr /= 0) return
            my_param1 = first_my_param1
         end if

         if (vary_my_param2) then
            call star_astero_procs% will_set_my_param( &
               s% id, 2, next_my_param2_to_try, ierr)
            if (ierr /= 0) return
            my_param2 = next_my_param2_to_try
         else
            call star_astero_procs% will_set_my_param( &
               s% id, 2, first_my_param2, ierr)
            if (ierr /= 0) return
            my_param2 = first_my_param2
         end if

         if (vary_my_param3) then
            call star_astero_procs% will_set_my_param( &
               s% id, 3, next_my_param3_to_try, ierr)
            if (ierr /= 0) return
            my_param3 = next_my_param3_to_try
         else
            call star_astero_procs% will_set_my_param( &
               s% id, 3, first_my_param3, ierr)
            if (ierr /= 0) return
            my_param3 = first_my_param3
         end if
         
         s% job% relax_initial_mass = .true.
         s% initial_mass = s% job% new_mass
         
         initial_Y = Y
         !s% initial_Z = Z << don't do this. it interferes with use of zams file.
         
         s% job% initial_h1 = X
         s% job% initial_h2 = 0
         s% job% initial_he3 = Y_frac_he3*Y
         s% job% initial_he4 = Y - s% job% initial_he3
         s% job% set_uniform_initial_composition = .true. 
         
         current_Y = Y
         current_FeH = FeH
         current_mass = s% job% new_mass
         current_alpha = s% mixing_length_alpha
         current_f_ov = f_ov

         current_my_param1 = my_param1
         current_my_param2 = my_param2
         current_my_param3 = my_param3
         
         current_h1 = X
         current_he3 = s% job% initial_he3
         current_he4 = s% job% initial_he4
         current_Z = Z

         if (f_ov > 0._dp) then
            s% overshoot_scheme(1) = 'exponential'
            s% overshoot_zone_type(1) = 'any'
            s% overshoot_zone_loc(1) = 'any'
            s% overshoot_bdy_loc(1) = 'any'
            s% overshoot_f(1) = f_ov
            s% overshoot_f0(1) = f0_ov_div_f_ov*f_ov
         else
            s% overshoot_scheme(1) = ''
            s% overshoot_zone_type(1) = ''
            s% overshoot_zone_loc(1) = ''
            s% overshoot_bdy_loc(1) = ''
            s% overshoot_f(1) = 0d0
            s% overshoot_f0(1) = 0d0
         end if
         
      end subroutine astero_extras_controls
      
      
      subroutine astero_extras_startup(id, restart, ierr)
         integer, intent(in) :: id
         logical, intent(in) :: restart
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         call star_astero_procs% extras_startup(id, restart, ierr)
         if (.not. restart) then
            call alloc_extra_info(s)
         else ! it is a restart
            call unpack_extra_info(s)
         end if
      end subroutine astero_extras_startup


      integer function astero_how_many_extra_history_columns(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         astero_how_many_extra_history_columns = &
            star_astero_procs% how_many_extra_history_columns(id)
         if (.not. s% job% astero_just_call_my_extras_check_model) &
            astero_how_many_extra_history_columns = &
               astero_how_many_extra_history_columns + num_extra_history_columns
      end function astero_how_many_extra_history_columns
      
      
      subroutine astero_data_for_extra_history_columns(id, n, astero_names, vals, ierr)
         integer, intent(in) :: id, n
         character (len=maxlen_history_column_name) :: astero_names(n)
         real(dp) :: vals(n)
         integer, intent(out) :: ierr
         integer :: i, num_extra
         type (star_info), pointer :: s
         
         include 'formats'
         
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         
         call star_astero_procs% data_for_extra_history_columns( &
            id, n, astero_names, vals, ierr)
         if (ierr /= 0) return
         if (s% job% astero_just_call_my_extras_check_model) return
         
         num_extra = star_astero_procs% how_many_extra_history_columns(id)
         
         i = num_extra+1
         astero_names(i) = 'chi2'         
         i = i+1
         astero_names(i) = 'delta_nu'         
         i = i+1
         astero_names(i) = 'delta_nu_model'         
         i = i+1
         astero_names(i) = trim(surf_coef1_name)
         i = i+1
         astero_names(i) = trim(surf_coef2_name)
         
         if (i /= (num_extra_history_columns + num_extra)) then
            write(*,2) 'i', i
            write(*,2) 'num_extra_history_columns', num_extra_history_columns
            call mesa_error(__FILE__,__LINE__,'bad num_extra_history_columns')
         end if
         
         i = num_extra+1
         vals(i) = chi2
         i = i+1
         vals(i) = delta_nu
         i = i+1
         vals(i) = delta_nu_model
         i = i+1
         vals(i) = surf_coef1
         i = i+1
         vals(i) = surf_coef2
         
         
      end subroutine astero_data_for_extra_history_columns

      
      integer function astero_how_many_extra_profile_columns(id)
         integer, intent(in) :: id
         astero_how_many_extra_profile_columns = &
            star_astero_procs% how_many_extra_profile_columns(id)
      end function astero_how_many_extra_profile_columns
      
      
      subroutine astero_data_for_extra_profile_columns( &
            id, n, nz, astero_names, vals, ierr)
         integer, intent(in) :: id, n, nz
         character (len=maxlen_profile_column_name) :: astero_names(n)
         real(dp) :: vals(nz,n)
         integer, intent(out) :: ierr
         integer :: k
         ierr = 0
         call star_astero_procs% data_for_extra_profile_columns( &
            id, n, nz, astero_names, vals, ierr)
      end subroutine astero_data_for_extra_profile_columns

      
      subroutine astero_extras_after_evolve(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         integer :: iounit, ckm
         type (star_info), pointer :: s
         character (len=256) :: filename
         include 'formats'
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         if ((eval_chi2_at_target_age_only .and. s% star_age >= age_target) .or. &
             (include_age_in_chi2_spectro .and. s% star_age >= max_age_for_chi2) .or. &
             save_info_for_last_model) then
            !write(*,*) 'call do_astero_extras_check_model before terminate'
            ckm = do_astero_extras_check_model(s, id)
            !write(*,*) 'done do_astero_extras_check_model before terminate'
         end if         
         call star_astero_procs% extras_after_evolve(id, ierr)
         if (save_info_for_last_model) then
            write(*,1) 'chi2', chi2
            call get_all_el_info(s,ierr)
            if (ierr /= 0) return
            call store_best_info(s)

            filename = trim(astero_results_directory) // '/' // trim(last_model_save_info_filename)

            iounit = alloc_iounit(ierr)
            if (ierr /= 0) return
            open(unit=iounit, file=filename, &
               action='write', status='replace', iostat=ierr)
            if (ierr /= 0) then
               write(*,'(a)') 'failed to open last_model_save_info_filename ' // filename
               return
            end if
            write(*,*) 'write ' // filename
            write(*,*) 'call show_best'
            call show_best(iounit)
            write(*,*) 'done show_best'
            call free_iounit(iounit)
         end if
         
         if (s% job% astero_just_call_my_extras_check_model) return

      end subroutine astero_extras_after_evolve
      
      
      ! routines for saving and restoring extra data so can do restarts
      
      subroutine alloc_extra_info(s)
         integer, parameter :: extra_info_alloc = 1
         type (star_info), pointer :: s
         call move_extra_info(s,extra_info_alloc)
      end subroutine alloc_extra_info
      
      
      subroutine unpack_extra_info(s)
         integer, parameter :: extra_info_get = 2
         type (star_info), pointer :: s
         call move_extra_info(s,extra_info_get)
      end subroutine unpack_extra_info
      
      
      subroutine store_extra_info(s)
         integer, parameter :: extra_info_put = 3
         type (star_info), pointer :: s
         call move_extra_info(s,extra_info_put)
      end subroutine store_extra_info
      
      
      subroutine move_extra_info(s,op)
         integer, parameter :: extra_info_alloc = 1
         integer, parameter :: extra_info_get = 2
         integer, parameter :: extra_info_put = 3
         type (star_info), pointer :: s
         integer, intent(in) :: op
         
         integer :: i, j, num_ints, num_dbls, ierr
         
         i = 0
         ! call move_int or move_flg 
         num_ints = i
         
         i = 0
         ! call move_dbl
         num_dbls = i
         
         if (op /= extra_info_alloc) return
         if (num_ints == 0 .and. num_dbls == 0) return
         
         ierr = 0
         call star_alloc_extras(s% id, num_ints, num_dbls, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in star_alloc_extras'
            write(*,*) 'alloc_extras num_ints', num_ints
            write(*,*) 'alloc_extras num_dbls', num_dbls
            call mesa_error(__FILE__,__LINE__)
         end if
         
         contains
         
         subroutine move_dbl(dbl)
            real(dp) :: dbl
            i = i+1
            select case (op)
            case (extra_info_get)
               dbl = s% extra_work(i)
            case (extra_info_put)
               s% extra_work(i) = dbl
            end select
         end subroutine move_dbl
         
         subroutine move_int(int)
            integer :: int
            i = i+1
            select case (op)
            case (extra_info_get)
               int = s% extra_iwork(i)
            case (extra_info_put)
               s% extra_iwork(i) = int
            end select
         end subroutine move_int
         
         subroutine move_flg(flg)
            logical :: flg
            i = i+1
            select case (op)
            case (extra_info_get)
               flg = (s% extra_iwork(i) /= 0)
            case (extra_info_put)
               if (flg) then
                  s% extra_iwork(i) = 1
               else
                  s% extra_iwork(i) = 0
               end if
            end select
         end subroutine move_flg
      
      end subroutine move_extra_info



      end module extras_support
