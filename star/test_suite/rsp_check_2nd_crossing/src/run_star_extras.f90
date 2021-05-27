! ***********************************************************************
!
!   Copyright (C) 2018-2019  The MESA Team
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

      implicit none
      
      include "test_suite_extras_def.inc"

      contains

      include "test_suite_extras.inc"
      
      
      subroutine extras_controls(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         s% extras_startup => extras_startup
         s% extras_check_model => extras_check_model
         s% extras_finish_step => extras_finish_step
         s% extras_after_evolve => extras_after_evolve
         s% how_many_extra_history_columns => how_many_extra_history_columns
         s% data_for_extra_history_columns => data_for_extra_history_columns
         s% how_many_extra_profile_columns => how_many_extra_profile_columns
         s% data_for_extra_profile_columns => data_for_extra_profile_columns  
         
         if (s% use_other_rsp_build_model) &
            s% other_rsp_build_model => rsp_check_2nd_crossing

      end subroutine extras_controls
      

      subroutine rsp_check_2nd_crossing(id, ierr)
         use interp_1d_def, only: pm_work_size
         use interp_1d_lib
         integer, intent(in) :: id
         integer, intent(out) :: ierr  
         type (star_info), pointer :: s
         
         integer, parameter :: io_in=34, io_out=35, max_n = 200, max_cnt = 9000
         real(dp) :: logT1, logL1, logT2, logL2, logT3, logL3, logT4, logL4, &
            delta_Teff, mass, X, Z, log_Teff, log_L, prev_Teff, &
            max_T, min_T, Teff_red_edge, Teff_blue_edge, offset
         integer :: col_model_number, col_star_age, col_log_Teff, col_log_L, &
            skip_cols, num_cols_to_read, i, n, cnt, num_beyond_blue_edge, &
            model_cnt, num_models, i_red, i_blue
         integer, allocatable, dimension(:) :: modnums, model
         real(dp), allocatable, dimension(:) :: &
            vals, growth, period, temp, lum, Ts, Ls, ages
         real(dp), pointer :: f(:,:)
         real(dp), pointer, dimension(:) :: &
            f1, work1, x_old, v_old, x_new, v_new
         logical :: okay, finished_1st_crossing, have_first, in_2nd_crossing, just_failed
         
         include 'formats'
         
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         write(*,*) 'rsp_check_2nd_crossing'
         
         delta_Teff = s% x_ctrl(1)
         ! approx edges
         logT1 = s% x_ctrl(2)
         logL1 = s% x_ctrl(3)
         logT2 = s% x_ctrl(4)
         logL2 = s% x_ctrl(5)
         logT3 = s% x_ctrl(6)
         logL3 = s% x_ctrl(7)
         logT4 = s% x_ctrl(8)
         logL4 = s% x_ctrl(9)
         
         skip_cols = s% x_integer_ctrl(1)
         col_model_number = s% x_integer_ctrl(2)
         col_star_age = s% x_integer_ctrl(3)
         col_log_Teff = s% x_integer_ctrl(4)
         col_log_L = s% x_integer_ctrl(5)
         
         num_cols_to_read = max(col_model_number, col_star_age, col_log_Teff, col_log_L)
         allocate(vals(num_cols_to_read), &
            growth(max_cnt), period(max_cnt), temp(max_cnt), lum(max_cnt), &
            Ts(max_cnt), Ls(max_cnt), ages(max_cnt), modnums(max_cnt), model(max_cnt), &
            f1(4*max_cnt), work1(max_cnt*pm_work_size), &
            x_old(max_n), v_old(max_n), x_new(max_n), v_new(max_n))
         
         open(unit=io_in, file=trim(s% x_character_ctrl(1)), status='old', action='read', iostat=ierr)
         if (ierr /= 0) then
            write(*,*) 'failed to open history data file ' // trim(s% x_character_ctrl(1))
            return
         end if        
         
         do i=1,skip_cols
            read(io_in,*)
         end do
         
         model_cnt = 0
         read_loop: do
            read(io_in,fmt=*,iostat=ierr) vals(1:num_cols_to_read)
            if (ierr /= 0) exit read_loop
            if (model_cnt >= max_cnt) stop 'need to increase max_cnt'
            model_cnt = model_cnt + 1
            modnums(model_cnt) = int(vals(col_model_number))
            ages(model_cnt) = vals(col_star_age)
            Ts(model_cnt) = exp10(vals(col_log_Teff))
            Ls(model_cnt) = exp10(vals(col_log_L))         
         end do read_loop
         num_models = model_cnt
         close(io_in)
         
         open(unit=io_out, file=TRIM(s% x_character_ctrl(2)), status='REPLACE', iostat=ierr)
         if (ierr /= 0) then
            write(*,*) 'failed to open output data file ' // trim(s% x_character_ctrl(2))
            return
         end if

         mass = s% RSP_mass
         X = s% RSP_X
         Z = s% RSP_Z
         
         write(io_out,1) 'RSP_mass', s% RSP_mass
         write(io_out,1) 'RSP_X', s% RSP_X
         write(io_out,1) 'RSP_Z', s% RSP_Z
         write(io_out,1) 'RSP_alfam', s% RSP_alfam
         write(io_out,1) 'RSP_alfap', s% RSP_alfap
         write(io_out,1) 'RSP_alfat', s% RSP_alfat
         write(io_out,1) 'RSP_gammar', s% RSP_gammar
         write(io_out,*)
         
         write(io_out,'(99a20)') 'model_number', 'period(d)', 'growth', &
            'Teff', 'L', 'star_age'

         s% RSP_nmodes = 1 ! just F
         finished_1st_crossing = .false.
         in_2nd_crossing = .false.
         have_first = .false.
         just_failed = .false.
         n = 0
         cnt = 0
         num_beyond_blue_edge = 0
         prev_Teff = 1d99
         search_loop: do model_cnt=1,num_models
            s% RSP_Teff = Ts(model_cnt)
            s% RSP_L = Ls(model_cnt)
            log_L = log10(s% RSP_L)
            if (.not. finished_1st_crossing) then
               min_T = exp10(get_red_logT(log_L))
               finished_1st_crossing = (s% RSP_Teff < min_T)
               prev_Teff = s% RSP_Teff
               cycle search_loop
            end if
            if (s% RSP_Teff < prev_Teff) then
               !write(*,2) 'still going to lower T', modnums(model_cnt), s% RSP_Teff, prev_Teff
               prev_Teff = s% RSP_Teff
               cycle search_loop ! still going to lower Ts
            end if
            max_T = exp10(get_blue_logT(log_L))            
            min_T = exp10(get_red_logT(log_L))            
            if (s% RSP_Teff < min_T - 4*delta_Teff) then
               !write(*,2) 'too far from red edge', modnums(model_cnt), s% RSP_Teff, min_T
               prev_Teff = s% RSP_Teff
               cycle search_loop ! too far from red edge
            end if
            if (have_first .and. s% RSP_Teff - delta_Teff < prev_Teff) then
               !write(*,2) 'too close to prev', modnums(model_cnt), s% RSP_Teff - delta_Teff - prev_Teff
               cycle search_loop ! too close to prev
            end if
            !if (s% RSP_Teff < 0.5d0*(min_T + max_T)) cycle search_loop ! too far from blue edge
            write(*,*) 'call star_do1_rsp_build model Teff L', &
               modnums(model_cnt), s% RSP_Teff, s% RSP_L
            call star_do1_rsp_build(s,ierr)
            if (ierr /= 0) then
               write(*,*) 'failed in star_do1_rsp_build'
               ierr = 0
               if (num_beyond_blue_edge > 0 .or. just_failed) exit search_loop
               just_failed = .true.
               cycle search_loop
               !stop 'failed in star_do1_rsp_build'
            end if
            just_failed = .false.
            n = n+1
            period(n) = s% rsp_LINA_periods(1)
            growth(n) = s% rsp_LINA_growth_rates(1)
            temp(n) = s% RSP_Teff
            lum(n) = s% RSP_L
            if ((.not. have_first) .and. growth(n) > 0d0) then
               write(*,*) 'failed to find red edge'
               ierr = -1
               return
            end if
            have_first = .true.
            model(n) = model_cnt
            if (growth(n) > 0d0) in_2nd_crossing = .true.
            if (growth(n) <= 0d0 .and. in_2nd_crossing) then
               num_beyond_blue_edge = num_beyond_blue_edge + 1
            end if
            write(io_out,'(i20,99(1pd20.10))') modnums(model_cnt), &
               period(n)/86400.d0, growth(n), s% RSP_Teff, s% RSP_L, ages(model_cnt)
            if (num_beyond_blue_edge == 2 .or. n == max_n) exit search_loop
            prev_Teff = s% RSP_Teff
         end do search_loop
         
         if (n == 0 .or. growth(1) >= 0d0 .or. num_beyond_blue_edge < 1) then
            write(*,2) 'n', n
            write(*,2) 'num_beyond_blue_edge', num_beyond_blue_edge
            write(*,1) 'growth(1)', growth(1)
            ierr = -1
            return
         end if
         
         Teff_red_edge = -1
         if (growth(2) < 0d0) then
            do i=4,n
               if (growth(i) > 0d0 .and. growth(i-1) > 0d0 .and. &
                   growth(i-2) < 0d0 .and. growth(i-3) < 0d0) then
                  i_red = i
                  x_old(1:4) = growth(i-3:i)
                  v_old(1:4) = temp(i-3:i)
                  x_new(1) = 0d0
                  !write(*,*) 'i_red', i_red
                  !write(*,*) 'x_old', x_old(1:4)
                  !write(*,*) 'v_old', v_old(1:4)
                  call interpolate_vector_pm( &
                     4, x_old, 1, x_new, v_old, v_new, work1, 'red edge', ierr)    
                  if (ierr /= 0) stop 'failed in interpolate_vector_pm red edge'
                  Teff_red_edge = v_new(1)
                  write(*,*)
                  write(*,'(a20,f20.10)') 'Teff_red_edge', Teff_red_edge
                  write(io_out,*)
                  write(io_out,'(a20,f20.10)') 'Teff_red_edge', Teff_red_edge
                  exit
               end if
            end do
         else if (growth(3) > 0d0) then
            i = 3
            i_red = i
            x_old(1:3) = growth(1:3)
            v_old(1:3) = temp(1:3)
            x_new(1) = 0d0
            !write(*,*) 'i_red', i_red
            !write(*,*) 'x_old', x_old(1:3)
            !write(*,*) 'v_old', v_old(1:3)
            call interpolate_vector_pm( &
               3, x_old, 1, x_new, v_old, v_new, work1, 'red edge', ierr)    
            if (ierr /= 0) stop 'failed in interpolate_vector_pm red edge'
            Teff_red_edge = v_new(1)
            write(*,*)
            write(*,'(a20,f20.10)') 'Teff_red_edge', Teff_red_edge
            write(io_out,*)
            write(io_out,'(a20,f20.10)') 'Teff_red_edge', Teff_red_edge
         end if
         if (Teff_red_edge < 0d0) stop 'failed to find red edge'
         
         Teff_blue_edge = -1
         if (num_beyond_blue_edge == 2) then
            do i=n-3,1,-1
               if (growth(i) > 0d0 .and. growth(i+1) > 0d0 .and. &
                   growth(i+2) < 0d0 .and. growth(i+3) < 0d0) then
                  i_blue = i
                  x_old(1:4) = growth(i:i+3)
                  v_old(1:4) = temp(i:i+3)
                  x_new(1) = 0d0
                  !write(*,*) 'i_blue', i_blue
                  !write(*,*) 'x_old', x_old(1:4)
                  !write(*,*) 'v_old', v_old(1:4)
                  call interpolate_vector_pm( &
                     4, x_old, 1, x_new, v_old, v_new, work1, 'blue edge', ierr)    
                  if (ierr /= 0) stop 'failed in interpolate_vector_pm blue edge'
                  Teff_blue_edge = v_new(1)
                  write(*,'(a20,f20.10)') 'Teff_blue_edge', Teff_blue_edge
                  write(io_out,'(a20,f20.10)') 'Teff_blue_edge', Teff_blue_edge
                  exit
               end if
            end do
         else if (num_beyond_blue_edge == 1 .and. n > 3) then
            if (growth(n) < 0d0 .and. growth(n-1) > 0d0 .and. growth(n-2) > 0d0) then
               i = n-2
               i_blue = i
               x_old(1:3) = growth(i:i+2)
               v_old(1:3) = temp(i:i+2)
               x_new(1) = 0d0
               !write(*,*) 'i_blue', i_blue
               !write(*,*) 'x_old', x_old(1:3)
               !write(*,*) 'v_old', v_old(1:3)
               call interpolate_vector_pm( &
                  3, x_old, 1, x_new, v_old, v_new, work1, 'blue edge', ierr)    
               if (ierr /= 0) stop 'failed in interpolate_vector_pm blue edge'
               Teff_blue_edge = v_new(1)
               write(*,'(a20,f20.10)') 'Teff_blue_edge', Teff_blue_edge
               write(io_out,'(a20,f20.10)') 'Teff_blue_edge', Teff_blue_edge
            end if
         end if
         write(*,*)
         write(io_out,*)
         if (Teff_blue_edge < 0d0) then
            write(*,*) 'failed to find blue edge'         
            return
         end if
         
         f(1:4,1:n) => f1(1:4*n)
         do i=1,n
            f(1,i) = lum(i)
         end do
         call interp_pm(temp, n, f1, pm_work_size, work1, 'rsp', ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in interp_pm Ts'
            return
         end if
         
         offset = 0d0
         write(*,'(a10,2a20)') 'offset', 'Teff', 'L'
         write(io_out,'(a10,2a20)') 'offset', 'Teff', 'L'
         do 
            x_new(1) = Teff_blue_edge - offset
            if (x_new(1) < Teff_red_edge) x_new(1) = Teff_red_edge
            call interp_values(temp, n, f1, 1, x_new, v_new, ierr)
            if (ierr /= 0) stop 'failed in interp_values Ts'
            write(*,'(i10,2f20.10)') int(Teff_blue_edge - x_new(1)), x_new(1), v_new(1)
            write(io_out,'(i10,2f20.10)') int(Teff_blue_edge - x_new(1)), x_new(1), v_new(1)
            if (x_new(1) == Teff_red_edge) exit
            offset = offset + 100d0
         end do
         write(io_out,*)
         write(*,*)
         
         close(io_out)

         write(*,*) 'done rsp_check_2nd_crossing'
         write(*,*) TRIM(s% x_character_ctrl(2))

         !deallocate(f1, work1, x_new, v_new)
         
         ierr = -1 ! to force termination of run
         
         contains
         
         real(dp) function get_blue_logT(log_L)
            real(dp), intent(in) :: log_L
            get_blue_logT = logT2 + (log_L - logL2)*(logT1 - logT2)/(logL1 - logL2)
         end function get_blue_logT
         
         real(dp) function get_red_logT(log_L)
            real(dp), intent(in) :: log_L
            get_red_logT = logT4 + (log_L - logL4)*(logT3 - logT4)/(logL3 - logL4)
         end function get_red_logT
         
      end subroutine rsp_check_2nd_crossing
      
      
      subroutine extras_startup(id, restart, ierr)
         integer, intent(in) :: id
         logical, intent(in) :: restart
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         call test_suite_startup(s, restart, ierr)
      end subroutine extras_startup
      

      ! returns either keep_going or terminate.
      integer function extras_finish_step(id)
         integer, intent(in) :: id
         integer :: ierr
         real(dp) :: target_period, rel_run_E_err
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_finish_step = keep_going
         if (s% x_integer_ctrl(1) <= 0) return
         if (s% rsp_num_periods < s% x_integer_ctrl(1)) return
         write(*,*)
         write(*,*)
         write(*,*)
         target_period = s% x_ctrl(1)
         rel_run_E_err = s% cumulative_energy_error/s% total_energy
         write(*,*) 'rel_run_E_err', rel_run_E_err
         if (s% total_energy /= 0d0 .and. abs(rel_run_E_err) > 1d-5) then
            write(*,*) '*** BAD rel_run_E_error ***', &
            s% cumulative_energy_error/s% total_energy
         else if (abs(s% rsp_period/(24*3600) - target_period) > 1d-2) then
            write(*,*) '*** BAD ***', s% rsp_period/(24*3600) - target_period, &
               s% rsp_period/(24*3600), target_period
         else
            write(*,*) 'good match for period', &
               s% rsp_period/(24*3600), target_period
         end if
         write(*,*)
         write(*,*)
         write(*,*)
         extras_finish_step = terminate
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
         how_many_extra_history_columns = 0
      end function how_many_extra_history_columns
      
      
      subroutine data_for_extra_history_columns(id, n, names, vals, ierr)
         integer, intent(in) :: id, n
         character (len=maxlen_history_column_name) :: names(n)
         real(dp) :: vals(n)
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
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
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
      end subroutine data_for_extra_profile_columns
      
      
            

      end module run_star_extras
      
