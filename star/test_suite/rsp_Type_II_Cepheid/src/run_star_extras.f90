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
      use auto_diff

      implicit none
      
      include "test_suite_extras_def.inc"
      logical :: need_to_write_LINA_data

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
         s% extras_start_step => extras_start_step
         s% extras_finish_step => extras_finish_step
         s% extras_after_evolve => extras_after_evolve
         s% how_many_extra_history_columns => how_many_extra_history_columns
         s% data_for_extra_history_columns => data_for_extra_history_columns
         s% how_many_extra_profile_columns => how_many_extra_profile_columns
         s% data_for_extra_profile_columns => data_for_extra_profile_columns  
      end subroutine extras_controls

      
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
            need_to_write_LINA_data = len_trim(s% x_character_ctrl(10)) > 0
         else ! it is a restart
            need_to_write_LINA_data = .false.
         end if
      end subroutine extras_startup
      

      integer function extras_start_step(id)
         integer, intent(in) :: id
         integer :: ierr, io, i
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_start_step = keep_going
         if (need_to_write_LINA_data) then
            io = 61
            open(io,file=trim(s% x_character_ctrl(10)),status='unknown')
            write(io, '(99d16.5)') s% RSP_mass, s% RSP_L, s% RSP_Teff, &
               (s% rsp_LINA_periods(i), s% rsp_LINA_growth_rates(i), i=1, s% RSP_nmodes)
            close(io)
            write(*,*) 'write ' // trim(s% x_character_ctrl(10))
            need_to_write_LINA_data = .false.
         end if
      end function extras_start_step
      

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
         write(*,'(A)')
         write(*,'(A)')
         write(*,'(A)')
         target_period = s% x_ctrl(1)
         rel_run_E_err = s% cumulative_energy_error/s% total_energy
         write(*,*) 'rel_run_E_err', rel_run_E_err
         if (s% total_energy /= 0d0 .and. abs(rel_run_E_err) > 1d-4) then
            write(*,*) '*** BAD rel_run_E_error ***', &
            s% cumulative_energy_error/s% total_energy
         else if (abs(s% rsp_period/(24*3600) - target_period) > 3d-2) then
            write(*,*) '*** BAD ***', s% rsp_period/(24*3600) - target_period, &
               s% rsp_period/(24*3600), target_period
         else
            write(*,*) 'good match for period', &
               s% rsp_period/(24*3600), target_period
         end if
         write(*,'(A)')
         write(*,'(A)')
         write(*,'(A)')
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
         names(i) = 'num_periods'; vals(i) = s% RSP_num_periods; i=i+1
         names(i) = 'period'; vals(i) = s% RSP_period/(24*3600); i=i+1
         names(i) = 'growth'; vals(i) = s% rsp_GREKM_avg_abs; i=i+1
         names(i) = 'max_v_div_cs'; vals(i) = 0; i=i+1
         names(i) = 'delta_R'; vals(i) = s% rsp_DeltaR; i=i+1
         names(i) = 'delta_Teff'; vals(i) = 0; i=i+1
         names(i) = 'delta_logL'; vals(i) = s% rsp_DeltaMag/2.5; i=i+1
         names(i) = 'delta_Mag'; vals(i) = 0; i=i+1
      end subroutine data_for_extra_history_columns

      
      integer function how_many_extra_profile_columns(id)
         use star_def, only: star_info
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_profile_columns = 1
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
         names(1) = 'sign_Lc'
         do k=1,nz
            if (abs(s% Lc(k)) < 1d-6) then
               vals(k,1) = 0d0
            else
               vals(k,1) = sign(1d0,s% Lc(k))
            end if
         end do
      end subroutine data_for_extra_profile_columns
      
      


      end module run_star_extras
      
