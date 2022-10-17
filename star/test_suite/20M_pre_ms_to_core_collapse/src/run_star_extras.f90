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
   use auto_diff
   use gyre_lib
   
   implicit none
   
   include "test_suite_extras_def.inc"
   
   !alpha_mlt_routine
   !alpha_H = s% x_ctrl(21)
   !alpha_other = s% x_ctrl(22)
   !H_limit = s% x_ctrl(23)

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
      s% other_alpha_mlt => alpha_mlt_routine
   end subroutine extras_controls
   
   
   subroutine alpha_mlt_routine(id, ierr)
      use chem_def, only : ih1
      integer, intent(in) :: id
      integer, intent(out) :: ierr
      type (star_info), pointer :: s
      integer :: k, h1
      real(dp) :: alpha_H, alpha_other, H_limit
      include 'formats'
      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return
      alpha_H = s% x_ctrl(21)
      alpha_other = s% x_ctrl(22)
      H_limit = s% x_ctrl(23)
      h1 = s% net_iso(ih1)
      !write(*,1) 'alpha_H', alpha_H
      !write(*,1) 'alpha_other', alpha_other
      !write(*,1) 'H_limit', H_limit
      !write(*,2) 'h1', h1
      !write(*,2) 's% nz', s% nz
      if (alpha_H <= 0 .or. alpha_other <= 0 .or. h1 <= 0) return
      do k = 1, s% nz
         if (s% xa(h1, k) >= H_limit) then
            s% alpha_mlt(k) = alpha_H
         else
            s% alpha_mlt(k) = alpha_other
         end if
         !write(*,2) 'alpha_mlt', k, s% alpha_mlt(k),
      end do
      !stop
   end subroutine alpha_mlt_routine
   
   
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
      real(dp) :: dt, m
      integer :: k, nz
      include 'formats'
      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return
      
      write(*, '(A)')
      select case (s% x_integer_ctrl(1))
      case (7)
         ! put target info in TestHub output
         testhub_extras_names(1) = 'fe_core_mass'; testhub_extras_vals(1) = s% fe_core_mass
         
         if(s% fe_core_mass < 1d0) then
            write(*, 1) "Bad fe_core_mass", s%fe_core_mass
         else
            if(s% fe_core_infall > s% fe_core_infall_limit) then
               write(*, '(a)') 'all values are within tolerance'
            else
               write(*, '(a)') "Bad fe core infall"
            end if
         end if
      end select
      
      call test_suite_after_evolve(s, ierr)
      if (ierr /= 0) return
   
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
      character (len = maxlen_history_column_name) :: names(n)
      real(dp) :: vals(n)
      integer, intent(out) :: ierr
      type (star_info), pointer :: s
      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return
   end subroutine data_for_extra_history_columns
   
   
   integer function how_many_extra_profile_columns(id)
      use star_def, only : star_info
      integer, intent(in) :: id
      integer :: ierr
      type (star_info), pointer :: s
      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return
      how_many_extra_profile_columns = 1
   end function how_many_extra_profile_columns
   
   
   subroutine data_for_extra_profile_columns(id, n, nz, names, vals, ierr)
      use star_def, only : star_info, maxlen_profile_column_name
      use const_def, only : dp
      integer, intent(in) :: id, n, nz
      character (len = maxlen_profile_column_name) :: names(n)
      real(dp) :: vals(nz, n)
      integer, intent(out) :: ierr
      type (star_info), pointer :: s
      integer :: k
      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return
      names(1) = 'zbar_div_abar'
      do k = 1, s% nz
         vals(k, 1) = s% zbar(k) / s% abar(k)
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
      
      if (extras_finish_step == terminate) &
         s% termination_code = t_extras_finish_step
   end function extras_finish_step


end module run_star_extras
      
