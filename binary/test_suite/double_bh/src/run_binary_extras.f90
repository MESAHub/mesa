! ***********************************************************************
!
!   Copyright (C) 2012  Bill Paxton
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
module run_binary_extras
   
   use star_lib
   use binary_lib
   use star_def
   use const_def
   use math_lib
   use binary_def
   
   implicit none
   
   include "binary_test_suite_extras_def.inc"

contains

include "binary_test_suite_extras.inc"
   
   subroutine extras_binary_controls(binary_id, ierr)
      integer :: binary_id
      integer, intent(out) :: ierr
      type (binary_info), pointer :: b
      ierr = 0
      
      call binary_ptr(binary_id, b, ierr)
      if (ierr /= 0) then
         write(*, *) 'failed in binary_ptr'
         return
      end if
      
      ! Set these function pointers to point to the functions you wish to use in
      ! your run_binary_extras. Any which are not set, default to a null_ version
      ! which does nothing.
      b% how_many_extra_binary_history_header_items => how_many_extra_binary_history_header_items
      b% data_for_extra_binary_history_header_items => data_for_extra_binary_history_header_items
      b% how_many_extra_binary_history_columns => how_many_extra_binary_history_columns
      b% data_for_extra_binary_history_columns => data_for_extra_binary_history_columns
      
      b% extras_binary_startup => extras_binary_startup
      b% extras_binary_start_step => extras_binary_start_step
      b% extras_binary_check_model => extras_binary_check_model
      b% extras_binary_finish_step => extras_binary_finish_step
      b% extras_binary_after_evolve => extras_binary_after_evolve
      
      ! Once you have set the function pointers you want, then uncomment this (or set it in your star_job inlist)
      ! to disable the printed warning message,
      b% warn_binary_extra = .false.
   
   end subroutine extras_binary_controls
   
   
   integer function how_many_extra_binary_history_header_items(binary_id)
      use binary_def, only : binary_info
      integer, intent(in) :: binary_id
      how_many_extra_binary_history_header_items = 0
   end function how_many_extra_binary_history_header_items
   
   subroutine data_for_extra_binary_history_header_items(&
      binary_id, n, names, vals, ierr)
      type (binary_info), pointer :: b
      integer, intent(in) :: binary_id, n
      character (len = maxlen_binary_history_column_name) :: names(n)
      real(dp) :: vals(n)
      integer, intent(out) :: ierr
      ierr = 0
      call binary_ptr(binary_id, b, ierr)
      if (ierr /= 0) then
         write(*, *) 'failed in binary_ptr'
         return
      end if
   end subroutine data_for_extra_binary_history_header_items
   
   
   integer function how_many_extra_binary_history_columns(binary_id)
      use binary_def, only : binary_info
      integer, intent(in) :: binary_id
      how_many_extra_binary_history_columns = 0
   end function how_many_extra_binary_history_columns
   
   subroutine data_for_extra_binary_history_columns(binary_id, n, names, vals, ierr)
      use const_def, only : dp
      type (binary_info), pointer :: b
      integer, intent(in) :: binary_id
      integer, intent(in) :: n
      character (len = maxlen_binary_history_column_name) :: names(n)
      real(dp) :: vals(n)
      integer, intent(out) :: ierr
      real(dp) :: beta
      ierr = 0
      call binary_ptr(binary_id, b, ierr)
      if (ierr /= 0) then
         write(*, *) 'failed in binary_ptr'
         return
      end if
   
   end subroutine data_for_extra_binary_history_columns
   
   
   integer function extras_binary_startup(binary_id, restart, ierr)
      type (binary_info), pointer :: b
      integer, intent(in) :: binary_id
      integer, intent(out) :: ierr
      logical, intent(in) :: restart
      call binary_ptr(binary_id, b, ierr)
      if (ierr /= 0) then ! failure in  binary_ptr
         return
      end if
      
      call test_suite_startup(b, restart, ierr)
      
      extras_binary_startup = keep_going
   end function extras_binary_startup
   
   integer function extras_binary_start_step(binary_id, ierr)
      type (binary_info), pointer :: b
      integer, intent(in) :: binary_id
      integer, intent(out) :: ierr
      
      extras_binary_start_step = keep_going
      call binary_ptr(binary_id, b, ierr)
      if (ierr /= 0) then ! failure in  binary_ptr
         return
      end if
   
   end function  extras_binary_start_step
   
   !Return either keep_going, retry or terminate
   integer function extras_binary_check_model(binary_id)
      type (binary_info), pointer :: b
      integer, intent(in) :: binary_id
      integer :: ierr
      call binary_ptr(binary_id, b, ierr)
      if (ierr /= 0) then ! failure in  binary_ptr
         return
      end if
      extras_binary_check_model = keep_going
   
   end function extras_binary_check_model
   
   
   ! returns either keep_going or terminate.
   ! note: cannot request retry; extras_check_model can do that.
   integer function extras_binary_finish_step(binary_id)
      type (binary_info), pointer :: b
      integer, intent(in) :: binary_id
      integer :: ierr
      call binary_ptr(binary_id, b, ierr)
      
      if (ierr /= 0) then ! failure in  binary_ptr
         return
      end if
      
      if (b% ignore_rlof_flag .and. &
         abs(log10(abs(b% s1% L_nuc_burn_total * Lsun / b% s1% L(1)))) < 0.005 .and. &
         b% s1% star_age > 1d2) then
         ! if here, primary reached thermal equilibrium (reached ZAMS), so activate RLOF
         ! this is the amount of overflow of a q=1 system at L2, anything more than this
         ! is too much
         if (b% rl_relative_gap(1) > 0.320819224) then
            extras_binary_finish_step = terminate
            write(*, *) "Terminate due to overflow of L2 at ZAMS"
            return
         else if (b% rl_relative_gap(1) > 0d0) then
            write(*, *) "model is overflowing at ZAMS"
            ! initially overflowing system will evolve rapidly to q=1, limit mdot
            ! for contact system to work and remove hard limit on mdot change
            b% max_implicit_abs_mdot = 1d-3
            b% s1% delta_mdot_hard_limit = -1
            b% s2% delta_mdot_hard_limit = -1
         else
            write(*, *) "model is not overflowing at ZAMS"
         end if
         call binary_set_ignore_rlof_flag(b% binary_id, .false., ierr)
         if (ierr /= 0) then
            return
         end if
         b% terminate_if_L2_overflow = .true.
         write(*, *) "Engage RLOF!"
      else if (b% ignore_rlof_flag .and. &
         (abs(log10(abs(b% s1% L_nuc_burn_total * Lsun / b% s1% L(1)))) > 0.005 .or. &
            b% s1% star_age < 1d2)) then
         ! if here, still not in ZAMS, keep period fixed
         call binary_set_period_eccentricity(b% binary_id, &
            b% initial_period_in_days * secday, 0d0, ierr)
         if (ierr /= 0) then
            return
         end if
      end if
      
      ! check if stars are evolving homogeneously
      if (b% s1% center_h1 > 1d-3) then
         if (b% s1% center_he4 - b% s1% surface_he4 > 0.2) then
            extras_binary_finish_step = terminate
            write(*, *) "Terminate due to primary not evolving homogeneously"
            return
         end if
      end if
      if (b% s2% center_h1 > 1d-3) then
         if (b% s2% center_he4 - b% s2% surface_he4 > 0.2) then
            extras_binary_finish_step = terminate
            write(*, *) "Terminate due to secondary not evolving homogeneously"
            return
         end if
      end if
      
      ! turn a star into point mass once it depletes helium
      ! terminate once both deplete helium
      if (b% point_mass_i == 0) then
         ! if modelling q=1 system then just terminate at He depletion
         if (b% s1% center_he4 < 1d-5) then
            call binary_set_point_mass_i(b% binary_id, 1, ierr)
            b% eq_initial_bh_mass = b% m(1)
            if (ierr /= 0) then
               return
            end if
         else if (b% s2% center_he4 < 1d-5) then
            call binary_set_point_mass_i(b% binary_id, 2, ierr)
            b% eq_initial_bh_mass = b% m(2)
            if (ierr /= 0) then
               return
            end if
         end if
      end if
      ! terminate when the other star also depletes helium
      if (b% point_mass_i == 1) then
         if (b% s2% center_he4 < 1d-5) then
            extras_binary_finish_step = terminate
            write(*, *) "Terminate due to helium depletion"
            return
         end if
      else if (b% point_mass_i == 2) then
         if (b% s1% center_he4 < 1d-5) then
            extras_binary_finish_step = terminate
            write(*, *) "Terminate due to helium depletion"
            return
         end if
      end if
      
      extras_binary_finish_step = keep_going
   
   end function extras_binary_finish_step
   
   subroutine extras_binary_after_evolve(binary_id, ierr)
      type (binary_info), pointer :: b
      integer, intent(in) :: binary_id
      integer, intent(out) :: ierr
      call binary_ptr(binary_id, b, ierr)
      if (ierr /= 0) then ! failure in  binary_ptr
         return
      end if
      
      call test_suite_after_evolve(b, ierr)
   
   end subroutine extras_binary_after_evolve

end module run_binary_extras
