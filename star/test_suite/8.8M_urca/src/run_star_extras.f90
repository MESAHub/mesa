! ***********************************************************************
!
!   Copyright (C) 2011  Bill Paxton
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


      ! these routines are called by the standard run_star check_model
      contains

      include "test_suite_extras.inc"

      real(dp) function center_avg_x(s,j)
         type (star_info), pointer :: s
         integer, intent(in) :: j
         real(dp) :: sum_x, sum_dq, dx, dq
         integer :: k
         sum_x = 0
         sum_dq = 0
         do k = s% nz, 1, -1
            dq = s% dq(k)
            dx = s% xa(j,k)*dq
            if (sum_dq+dq >= s% center_avg_value_dq) then
               sum_x = sum_x+ dx*(s% center_avg_value_dq - sum_dq)/dq
               sum_dq = s% center_avg_value_dq
               exit
            end if
            sum_x = sum_x + dx
            sum_dq = sum_dq + dq
         end do
         center_avg_x = sum_x/sum_dq
      end function center_avg_x


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
      end subroutine extras_startup


      subroutine extras_after_evolve(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         real(dp) :: dt
         character (len=strlen) :: test
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         call test_suite_after_evolve(s, ierr)
      end subroutine extras_after_evolve


      ! returns either keep_going, retry, or terminate.
      integer function extras_check_model(id)

         use chem_def, only: img25, ina25, ina23, ine23

         integer, intent(in) :: id
         
         integer :: mg25, na25, na23, ne23, ierr
         real(dp) :: center_mg25, center_na25, center_na23, center_ne23
         type (star_info), pointer :: s

         ! test parameters
         real(dp), parameter :: x_limit = 0.01d0
         real(dp), parameter :: log_center_density_limit = 9.29d0
         real(dp), parameter :: log_center_temperature_limit = 8.75d0

         ! for debugging
         logical, parameter :: dbg = .false.
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         ! by default, keep going
         extras_check_model = keep_going

         ! get isotopes
         mg25 = s% net_iso(img25)
         na25 = s% net_iso(ina25)
         na23 = s% net_iso(ina23)
         ne23 = s% net_iso(ine23)

         center_mg25 = center_avg_x(s,mg25)
         center_na25 = center_avg_x(s,na25)
         center_na23 = center_avg_x(s,na23)
         center_ne23 = center_avg_x(s,ne23)

         if (dbg) write(*,'(6(F8.3))') &
              center_mg25, center_na25, center_na23, center_ne23, &
              s% log_center_temperature, s% log_center_density

         if (s% log_center_density >= log_center_density_limit) then
            extras_check_model = terminate

            if (center_mg25 > x_limit) then
               termination_code_str(t_xtra1) = 'FAIL: electron captures on Mg-25 failed to occur'
               s% termination_code = t_xtra1
               return
            endif

            if (center_na23 > x_limit) then
               termination_code_str(t_xtra2) = 'FAIL: electron captures on Na-23 failed to occur'
               s% termination_code = t_xtra2
               return
            endif

            if (s% log_center_temperature > log_center_temperature_limit ) then
               termination_code_str(t_xtra3) = 'FAIL: Urca process failed to produce expected cooling'
               s% termination_code = t_xtra3
               return
            else
               termination_code_str(t_xtra4) = 'PASS: Urca process cooling occurred as expected'
               s% termination_code = t_xtra4
               return
            endif
         endif

         return

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
         names(1) = 'zbar_div_abar'
         do k=1,s% nz
            vals(k,1) = s% zbar(k)/s% abar(k)
         end do
      end subroutine data_for_extra_profile_columns


      ! returns either keep_going or terminate.
      integer function extras_finish_step(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_finish_step = keep_going
      end function extras_finish_step



      end module run_star_extras
