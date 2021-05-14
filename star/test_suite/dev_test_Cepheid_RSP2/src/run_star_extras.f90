! ***********************************************************************
!
!   Copyright (C) 2010  Bill Paxton
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
 
      module run_star_extras

      use star_lib
      use star_def
      use const_def
      use math_lib
      use run_star_support
      
      implicit none
      
      include 'test_suite_extras_def.inc'

      ! summary info at time of recently completely period
      integer :: num_periods, run_num_steps_end_prev, &
         run_num_iters_end_prev, run_num_retries_end_prev
      real(dp) :: period, KE_growth, KE_growth_avg, prev_KE_max, &
         period_max_vsurf_div_cs, period_delta_R, period_delta_Teff, &
         period_delta_logL, period_delta_Mag
      ! info for period in progress
      real(dp) :: time_started, v_div_cs_max, &
         KE_min, KE_max, R_min, R_max, L_min, L_max, T_min, T_max
            
      contains

      include 'test_suite_extras.inc'
      
      
      subroutine extras_controls(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         s% extras_startup => extras_startup
         s% extras_start_step => extras_start_step
         s% extras_check_model => extras_check_model
         s% extras_finish_step => extras_finish_step
         s% extras_after_evolve => extras_after_evolve
         s% how_many_extra_history_columns => how_many_extra_history_columns
         s% data_for_extra_history_columns => data_for_extra_history_columns
         s% how_many_extra_profile_columns => how_many_extra_profile_columns
         s% data_for_extra_profile_columns => data_for_extra_profile_columns  
         s% other_photo_write => photo_write
         s% other_photo_read => photo_read
      end subroutine extras_controls


      subroutine photo_write(id, iounit)
         integer, intent(in) :: id, iounit
         write(iounit) num_periods, run_num_steps_end_prev, &
            run_num_iters_end_prev, run_num_retries_end_prev, &
            period, KE_growth, KE_growth_avg, prev_KE_max, &
            period_max_vsurf_div_cs, period_delta_R, period_delta_Teff, &
            period_delta_logL, period_delta_Mag, time_started, v_div_cs_max, &
            KE_min, KE_max, R_min, R_max, L_min, L_max, T_min, T_max
      end subroutine photo_write


      subroutine photo_read(id, iounit, ierr)
         integer, intent(in) :: id, iounit
         integer, intent(out) :: ierr
         ierr = 0
         read(iounit, iostat=ierr) num_periods, run_num_steps_end_prev, &
            run_num_iters_end_prev, run_num_retries_end_prev, &
            period, KE_growth, KE_growth_avg, prev_KE_max, &
            period_max_vsurf_div_cs, period_delta_R, period_delta_Teff, &
            period_delta_logL, period_delta_Mag, time_started, v_div_cs_max, &
            KE_min, KE_max, R_min, R_max, L_min, L_max, T_min, T_max
      end subroutine photo_read
      
      
      subroutine extras_startup(id, restart, ierr)
         integer, intent(in) :: id
         logical, intent(in) :: restart
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         call test_suite_startup(s, restart, ierr)    
         if (.not. restart) then     
            num_periods = 0
            run_num_steps_end_prev = 0
            run_num_iters_end_prev = 0
            run_num_retries_end_prev = 0
            period = 0
            KE_growth = 0
            KE_growth_avg = 0
            prev_KE_max = 0
            period_max_vsurf_div_cs = 0
            period_delta_R = 0
            period_delta_Teff = 0
            period_delta_logL = 0
            period_delta_Mag = 0
            time_started = 0
            v_div_cs_max = 0
            KE_min = 0
            KE_max = 0
            R_min = 0
            R_max = 0
            L_min = 0
            L_max = 0
            T_min = 0
            T_max = 0
         end if
      end subroutine extras_startup


      integer function extras_start_step(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         include 'formats'
         ierr = 0
         extras_start_step = keep_going
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
      end function extras_start_step
      
      
      ! returns either keep_going, retry, or terminate.
      integer function extras_finish_step(id)
         use chem_def
         integer, intent(in) :: id
         type (star_info), pointer :: s
         integer :: ierr, test_period
         real(dp) :: target_period
         logical :: doing_pulses
         include 'formats'
         
         extras_finish_step = terminate
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_finish_step = keep_going
                  
         doing_pulses = s% x_logical_ctrl(7)
         if (.not. doing_pulses) return
         target_period = s% x_ctrl(7)
         if (target_period <= 0d0) return
         if (.not. get_period_info()) return
         
         test_period = s% x_integer_ctrl(7)
         if (num_periods < test_period .or. test_period <= 0) return
         
         ! have finished test run
         call report_test_results
         extras_finish_step = terminate
         
         contains
         
         logical function get_period_info()
            real(dp) :: v_surf, v_surf_start, KE, KE_avg, min_period, time_ended, &
               delta_R, min_deltaR_for_periods, KE_growth_avg_frac_new, &
               min_period_div_target, cs
            include 'formats'
            get_period_info = .false.
         
            if (s% r(1) < R_min) R_min = s% r(1)
            if (s% r(1) > R_max) R_max = s% r(1)
            if (s% L(1) < L_min) L_min = s% L(1)
            if (s% L(1) > L_max) L_max = s% L(1)
            if (s% Teff < T_min) T_min = s% Teff
            if (s% Teff > T_max) T_max = s% Teff
            KE = s% total_radial_kinetic_energy_end
            if (KE > KE_max) KE_max = KE
            if (KE < KE_min) KE_min = KE
            
            if (s% v_flag) then
               v_surf = s% v(1)
               v_surf_start = s% v_start(1)
            else if (s% u_flag) then
               v_surf = s% u_face_val(1)
               v_surf_start = s% u_face_start(1)
            else
               stop 'extras_finish_step: both v_flag and u_flag are false'
            end if
            cs = s% csound(1)
            if (v_surf > v_div_cs_max*cs) v_div_cs_max = v_surf/cs
               
            ! period is completed when v_surf goes from positive to negative during step
            if (v_surf > 0d0 .or. v_surf_start < 0d0) return
            
            if (time_started == 0) then ! start of 1st cycle
               time_started = s% time
               run_num_steps_end_prev = s% model_number
               run_num_iters_end_prev = s% total_num_solver_iterations
               run_num_retries_end_prev = s% num_retries
               prev_KE_max = 0d0
               call init_min_max_info
               write(*,*) 'first maximum radius, period calculations starting at model, day', &
                  s% model_number, s% time/(24*3600)
               return
            end if
            
            delta_R = R_max - R_min
            min_deltaR_for_periods = s% x_ctrl(8)*Rsun
            if (min_deltaR_for_periods > 0d0) then
               if (delta_R < min_deltaR_for_periods) return ! filter out glitches
            end if
         
            time_ended = s% time
            if (abs(v_surf - v_surf_start) > 1d-10) & ! tweak the end time to match when v_surf == 0
               time_ended = s% time - v_surf*s% dt/(v_surf - v_surf_start)
            min_period_div_target = s% x_ctrl(10)
            min_period = target_period*(24*3600)*min_period_div_target
            if (min_period > 0d0 .and. &
                time_ended - time_started < min_period) return ! filter out glitches

            period = time_ended - time_started
            num_periods = num_periods + 1
         
            if (num_periods > 1) then
               KE_avg = 0.5d0*(KE_max + prev_KE_max)
               KE_growth = (KE_max - prev_KE_max)/KE_avg
               KE_growth_avg_frac_new = s% x_ctrl(9)
               KE_growth_avg = KE_growth_avg_frac_new*KE_growth + &
                  (1d0 - KE_growth_avg_frac_new)*KE_growth_avg
            end if
         
            period_delta_Teff = T_max - T_min
            period_delta_R = R_max - R_min
            period_delta_logL = log10(L_max/L_min)
            period_delta_Mag = 2.5d0*period_delta_logL
            period_max_vsurf_div_cs = v_div_cs_max
            prev_KE_max = KE_max
            write(*,'(i4,a14,i4,2(a14,f8.3),99(a14,f12.5))')  &
               num_periods, &
               'steps/cycle', s% model_number - run_num_steps_end_prev, &
               'iters/step',  &
                  dble(s% total_num_solver_iterations - run_num_iters_end_prev)/ &
                  dble(s% model_number - run_num_steps_end_prev), &
               'period (d)', period/(24*3600), &
               'KE growth', KE_growth_avg, &
               'delta R/Rsun', period_delta_R/Rsun, &
               'delta logL', period_delta_logL, &
               'delta Teff', period_delta_Teff, &
               'max vsurf/cs', period_max_vsurf_div_cs

            time_started = time_ended
            run_num_steps_end_prev = s% model_number
            run_num_iters_end_prev = s% total_num_solver_iterations
            run_num_retries_end_prev = s% num_retries
            call init_min_max_info
            get_period_info = .true.

         end function get_period_info
         
         subroutine init_min_max_info
            v_div_cs_max = 0d0
            KE_min = 1d99
            KE_max  = -1d99
            R_min = 1d99
            R_max = -1d99
            L_min = 1d99
            L_max = -1d99
            T_min = 1d99
            T_max = -1d99
         end subroutine init_min_max_info
         
         subroutine report_test_results
            real(dp) :: rel_run_E_err
            write(*,*)
            write(*,*)
            write(*,*)
            rel_run_E_err = s% cumulative_energy_error/s% total_energy
            write(*,*) 'rel_run_E_err', rel_run_E_err
            if (s% total_energy /= 0d0 .and. abs(rel_run_E_err) > 1d-5) then
               write(*,*) '*** BAD rel_run_E_error ***', &
               s% cumulative_energy_error/s% total_energy
            else if (abs(period/(24*3600) - target_period) > 1d-2) then
               write(*,*) '*** BAD period ***', period/(24*3600) - target_period, &
                  period/(24*3600), target_period
            else
               write(*,*) 'good match for period', &
                  period/(24*3600), target_period
            end if
            write(*,*)
            write(*,*)
            write(*,*)
         end subroutine report_test_results
         
      end function extras_finish_step
      
      
      subroutine extras_after_evolve(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         real(dp) :: dt
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         call test_suite_after_evolve(s, ierr)
      end subroutine extras_after_evolve


      ! returns either keep_going, retry, or terminate.
      integer function extras_check_model(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_check_model = keep_going         
      end function extras_check_model


      integer function how_many_extra_history_columns(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_history_columns = 8
      end function how_many_extra_history_columns
      
      
      subroutine data_for_extra_history_columns(id, n, names, vals, ierr)
         integer, intent(in) :: id, n
         character (len=maxlen_history_column_name) :: names(n)
         real(dp) :: vals(n)
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer :: i
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         i = 1
         names(i) = 'num_periods'; vals(i) = num_periods; i=i+1
         names(i) = 'period'; vals(i) = period/(24*3600); i=i+1
         names(i) = 'growth'; vals(i) = KE_growth_avg; i=i+1
         names(i) = 'max_v_div_cs'; vals(i) = period_max_vsurf_div_cs; i=i+1
         names(i) = 'delta_R'; vals(i) = period_delta_R/Rsun; i=i+1
         names(i) = 'delta_Teff'; vals(i) = period_delta_Teff; i=i+1
         names(i) = 'delta_logL'; vals(i) = period_delta_logL; i=i+1
         names(i) = 'delta_Mag'; vals(i) = period_delta_Mag; i=i+1
      end subroutine data_for_extra_history_columns


      integer function how_many_extra_profile_columns(id)
         use star_def, only: star_info
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_profile_columns = 0
      end function how_many_extra_profile_columns
      
      
      subroutine data_for_extra_profile_columns(id, n, nz, names, vals, ierr)
         use star_def, only: star_info, maxlen_profile_column_name
         use const_def, only: dp
         integer, intent(in) :: id, n, nz
         character (len=maxlen_profile_column_name) :: names(n)
         real(dp) :: vals(nz,n)
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer :: k
         include 'formats'
         ierr = 0            
      end subroutine data_for_extra_profile_columns


      end module run_star_extras
      
