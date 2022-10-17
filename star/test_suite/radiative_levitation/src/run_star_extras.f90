! ***********************************************************************
!
!   Copyright (C) 2011  The MESA Team
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
   
   
   ! these routines are called by the standard run_star check_model
contains

include "test_suite_extras.inc"
   
   
   subroutine set_op_mono_factors(id, ierr)
      use chem_def, only : ife56
      integer, intent(in) :: id
      integer, intent(out) :: ierr
      type (star_info), pointer :: s
      integer :: j
      include 'formats'
      ierr = 0
      call get_star_ptr(id, s, ierr)
      if (ierr /= 0) return
      s% op_mono_factors(:) = 1
      j = s% net_iso(ife56)
      if (j == 0) then
         write(*, *) 'failed to find fe56 in current net'
         ierr = -1
         return
      end if
      s% op_mono_factors(j) = s% x_ctrl(1)
      if (abs(1d0 - s% x_ctrl(1)) > 1d-6) write(*, 1) &
         'set_op_mono_factors -- increase fe56 opacity by factor ', s% x_ctrl(1)
   end subroutine set_op_mono_factors
   
   
   subroutine extras_controls(id, ierr)
      integer, intent(in) :: id
      integer, intent(out) :: ierr
      type (star_info), pointer :: s
      logical :: dir_exists
      character (len = 500) :: fname
      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return
      
      if (s% x_ctrl(1) > 0d0) then
         if (len_trim(s% op_mono_data_path) == 0) then
            ierr = -1
         else
            fname = trim(s% op_mono_data_path) // '/a06.140'
            inquire(file = trim(fname), exist = dir_exists)
            if (.not. dir_exists) then
               write(*, '(a)') ' control op_mono_data_path = "' // trim(s% op_mono_data_path) // '"'
               write(*, '(a)') ' the file "' // trim(fname) // '" does not exist, so skip this test.'
               ierr = -1
            end if
         end if
         if (ierr /= 0) then
            write(*, *) 'this test was intentionally skipped'
            write(*, *) 'pretend found expected effects of radiative levitation.'
         end if
         write(*, *) 'extras_controls set pointer to set_op_mono_factors'
         s% set_op_mono_factors => set_op_mono_factors
      end if
      
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
      use chem_lib, only : chem_get_iso_id
      integer, intent(in) :: id
      integer, intent(out) :: ierr
      type (star_info), pointer :: s
      real(dp) :: dt
      integer :: k, nz, species, fe56, ni58
      real(dp) :: sum_dq
      include 'formats'
      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return
      
      if (s% x_ctrl(1) > 0d0) then
         if (s% star_age < s% max_age) then
            write(*, 1) 'failed to reach max age', s% star_age, s% max_age
         else
            nz = s% nz
            species = s% species
            sum_dq = 0d0
            do k = 1, nz
               if (sum_dq > 1d-10) exit
               sum_dq = sum_dq + s% dq(k)
            end do
            if (s% mixing_type(k) /= convective_mixing) then
               write(*, 4) 'bad mixing type', k, s% mixing_type(k), convective_mixing
               return
            end if
            fe56 = s% net_iso(chem_get_iso_id('fe56'))
            ni58 = s% net_iso(chem_get_iso_id('ni58'))
            if (fe56 <= 0 .or. ni58 <= 0) then
               write(*, 4) 'missing fe56 or ni58', k, fe56, ni58
               return
            end if
            if (s% xa(fe56, k) < 2d-2 .or. s% xa(ni58, k) < 7d-3) then
               write(*, 2) 'too little fe56 or ni58', k, s% xa(fe56, k), s% xa(ni58, k)
               return
            end if
            do k = 1, nz
               if ((s% xa(fe56, k) + s% xa(ni58, k)) > 0.8) then
                  write(*, 2) 'a region dominated by fe56/ni58 has formed', k, s% xa(fe56, k), s% xa(ni58, k)
                  return
               end if
            end do
            write(*, *) 'found expected effects of radiative levitation'
         end if
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
      how_many_extra_profile_columns = 0
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
   end subroutine data_for_extra_profile_columns
   
   
   ! returns either keep_going or terminate.
   integer function extras_finish_step(id)
      integer, intent(in) :: id
      integer :: ierr
      type (star_info), pointer :: s
      include 'formats'
      
      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return
      
      extras_finish_step = keep_going
   
   end function extras_finish_step


end module run_star_extras
      
