! ***********************************************************************
!
!   Copyright (C) 2010-2019  The MESA Team
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
      use utils_lib, only: mesa_error
            
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
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         real(dp) :: dt
         real(dp) :: max_abs_v
         include 'formats'
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         if (s% u_flag) then
            max_abs_v = maxval(abs(s% u(1:s% nz)/s% csound(1:s% nz)))
         else
            max_abs_v = maxval(abs(s% v(1:s% nz)/s% csound(1:s% nz)))
         end if
         if (max_abs_v > 0.1d0) then
            write(*,1) 'max_abs_v_div_cs is too large', max_abs_v
         else if (s% time < 0.99d5) then
            write(*,1) 'stop time should be 1d5 seconds', s% time
         else
            write(*,1) 'max_abs_v_div_cs is small enough', max_abs_v
         end if
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
         integer :: k
         real(dp) :: t
         type (star_info), pointer :: s
         include 'formats'
         ierr = 0
         if (n == 0) return
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         
         names(1) = 'Pmax_P'
         names(2) = 'Pmax_r_1m13'
         names(3) = 'Pmax_v'
         names(4) = 'Pmax_rho'
         names(5) = 'Pmax_T'
         names(6) = 'log_Pmax_P'
         names(7) = 'log_Pmax_r'
         names(8) = 'log_Pmax_v'
         names(9) = 'log_Pmax_rho'
         names(10) = 'log_Pmax_T'
         names(11) = 'Pmax_r_div_t'
         names(12) = 'Pmax_m_div_Msun'
         
         k = maxloc(s% Peos(1:s% nz), dim=1)
         if (k == s% nz) k = s% nz-1
         t = s% time
         
         !write(*,2) 'Pmax k r*1d-13', k, 0.5d0*(s% r(k)+s% r(k+1))*1d-13
                  
         vals(1) = s% Peos(k) ! Pmax_P
         vals(2) = 0.5d0*(s% r(k)+s% r(k+1))*1d-13 ! Pmax_r_1m13
         vals(3) = 0.5d0*(s% v(k)+s% v(k+1)) ! Pmax_v
         vals(4) = s% rho(k) ! Pmax_rho
         vals(5) = s% T(k) ! Pmax_T
         vals(6) = log10(s% Peos(k)) ! Pmax_P
         vals(7) = log10(0.5d0*(s% r(k)+s% r(k+1))) ! Pmax_r
         vals(8) = log10(0.5d0*(s% v(k)+s% v(k+1))) ! Pmax_v
         vals(9) = log10(s% rho(k)) ! Pmax_rho
         vals(10) = log10(s% T(k)) ! Pmax_T
         vals(11) = 0.5d0*(s% r(k)+s% r(k+1))/s% time ! Pmax_r_div_t
         vals(12) = s% m(k)/Msun ! Pmax_m_div_Msun
         
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
         use const_def, only: dp, avo, kerg
         integer, intent(in) :: id, n, nz
         character (len=maxlen_profile_column_name) :: names(n)
         real(dp) :: vals(nz,n)
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer :: k
         real(dp) :: gamma, energy
         ierr = 0
         if (n == 0) return
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
      end subroutine data_for_extra_profile_columns
      

      ! returns either keep_going or terminate.
      ! note: cannot request retry; extras_check_model can do that.
      integer function extras_finish_step(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         integer :: k
         include 'formats'
         extras_finish_step = keep_going
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         if (s% x_ctrl(9) <= 0) return
         do k = 1, s% nz
            if (s% v(k) > s% csound(k) .and. s% r(k) > s% x_ctrl(9)*Rsun) then
               write(*,1) 'shock has reached target location', s% r(k)/Rsun
               extras_finish_step = terminate
               return
            end if
         end do
      end function extras_finish_step
      
      


      end module run_star_extras
      
