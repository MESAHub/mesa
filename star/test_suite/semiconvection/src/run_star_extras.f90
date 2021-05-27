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
         use num_lib
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         real(dp) :: dt
         integer :: k, k1, k2, k3, k4
         type (star_info), pointer :: s
         logical :: okay
         
         include 'formats'

         
         okay = .true.

         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         call test_suite_after_evolve(s, ierr)
         
         write(*,*)
         k1 = 0
         k2 = 0
         k3 = 0
         k4 = 0
         do k = s% nz, 1, -1
            if (s% m(k) > 0.125d0*Msun .and. k1 == 0) k1 = k
            if (s% m(k) > 0.135d0*Msun .and. k2 == 0) k2 = k
            if (s% m(k) > 0.145d0*Msun .and. k3 == 0) k3 = k
         end do
         write(*,*)
         call check('mixing type at 0.125 Msun', dble(s% mixing_type(k1)), &
            dble(convective_mixing), dble(convective_mixing))
         call check('mixing type at 0.135 Msun', dble(s% mixing_type(k2)), &
            dble(semiconvective_mixing), dble(semiconvective_mixing))
         call check('mixing type at 0.145 Msun', dble(s% mixing_type(k3)), &
            dble(no_mixing), dble(no_mixing))
         call check('logT', avg_val(s% lnT)/ln10, 7.15d0, 7.31d0)
         call check('logRho', avg_val(s% lnd)/ln10, 1.75d0, 1.80d0)
         write(*,*)
         if (okay) write(*,'(a)') 'all values are within tolerances'
         write(*,*)
         
         
         contains
         
         real(dp) function avg_val(v)
            real(dp) :: v(:)
            avg_val = dot_product(v(k3:k2), s% dq(k3:k2)) / sum(s% dq(k3:k2))
         end function avg_val
         
         subroutine check(str, val, low, hi)
            real(dp), intent(in) :: val, low, hi
            character (len=*) :: str
            include 'formats'
            if (low <= val .and. val <= hi) then
               write(*,1) trim(str), val, low, hi
            else
               write(*,1) '*** BAD *** ' // trim(str), val, low, hi
               okay = .false.
            end if
         end subroutine check
         
         
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
      
