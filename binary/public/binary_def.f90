! ***********************************************************************
!
!   Copyright (C) 2013  Pablo Marchant
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
 

      module binary_def

      use star_lib
      use star_def
      use const_def
      
      implicit none

      real(dp) :: initial_binary_period ! (seconds)
      real(dp) :: min_binary_period ! (seconds)
      
      real(dp) :: initial_mass(2) ! (msun)

      integer, parameter :: maxlen_binary_history_column_name = 80
      integer, parameter :: binary_num_xtra_vals = 30
      
      ! time_step limit identifiers
      integer, parameter :: b_Tlim_comp = 1
      integer, parameter :: b_Tlim_roche = b_Tlim_comp + 1
      integer, parameter :: b_Tlim_jorb = b_Tlim_roche + 1
      integer, parameter :: b_Tlim_env = b_Tlim_jorb + 1
      integer, parameter :: b_Tlim_sep = b_Tlim_env + 1
      integer, parameter :: b_Tlim_ecc = b_Tlim_sep + 1
      integer, parameter :: b_Tlim_dm = b_Tlim_ecc + 1
      integer, parameter :: b_numTlim = b_Tlim_dm

      character (len=24) :: binary_dt_why_str(b_numTlim) ! indicates the reson for the timestep choice
         
      !interfaces for procedure pointers
      abstract interface

         subroutine other_rlo_mdot_interface(binary_id, rlo_mdot, ierr)
            use const_def, only: dp
            integer, intent(in) :: binary_id
            real(dp), intent(out) :: rlo_mdot
            integer, intent(out) :: ierr
         end subroutine other_rlo_mdot_interface
         
         integer function other_check_implicit_rlo_interface(binary_id, new_mdot)
            use const_def, only: dp
            integer, intent(in) :: binary_id
            real(dp), intent(out) :: new_mdot
         end function other_check_implicit_rlo_interface

         subroutine other_implicit_function_to_solve_interface(binary_id, function_to_solve, use_sum, detachment, ierr)
            use const_def, only: dp
            integer, intent(in) :: binary_id
            real(dp), intent(out) :: function_to_solve
            logical, intent(out) :: use_sum, detachment
            integer, intent(out) :: ierr
         end subroutine other_implicit_function_to_solve_interface

         subroutine other_tsync_interface(id, sync_type, Ftid, qratio, m, r_phot, osep, t_sync, ierr)
            use const_def, only: dp, strlen
            integer, intent(in) :: id
            character (len=strlen), intent(in) :: sync_type
            real(dp), intent(in) :: Ftid
            real(dp), intent(in) :: qratio
            real(dp), intent(in) :: m
            real(dp), intent(in) :: r_phot
            real(dp), intent(in) :: osep
            real(dp), intent(out) :: t_sync
            integer, intent(out) :: ierr
         end subroutine other_tsync_interface

         subroutine other_sync_spin_to_orbit_interface(id, nz, osep, qratio, rl, dt_next, Ftid, sync_type, sync_mode, ierr)
            use const_def, only: dp, strlen
            integer, intent(in) :: id
            integer, intent(in) :: nz
            real(dp), intent(in) :: osep
            real(dp), intent(in) :: qratio
            real(dp), intent(in) :: rl
            real(dp), intent(in) :: dt_next
            real(dp), intent(in) :: Ftid
            character (len=strlen), intent(in) :: sync_type
            character (len=strlen), intent(in) :: sync_mode
            integer, intent(out) :: ierr
         end subroutine other_sync_spin_to_orbit_interface

         subroutine other_mdot_edd_interface(binary_id, mdot_edd, mdot_edd_eta, ierr)
            use const_def, only: dp
            integer, intent(in) :: binary_id
            real(dp), intent(out) :: mdot_edd
            real(dp), intent(out) :: mdot_edd_eta
            integer, intent(out) :: ierr
         end subroutine other_mdot_edd_interface

         subroutine other_adjust_mdots_interface(binary_id, ierr)
            integer, intent(in) :: binary_id
            integer, intent(out) :: ierr
         end subroutine other_adjust_mdots_interface

         subroutine other_accreted_material_j_interface(binary_id, ierr)
            integer, intent(in) :: binary_id
            integer, intent(out) :: ierr
         end subroutine other_accreted_material_j_interface

         subroutine other_jdot_interface(binary_id, ierr)
            integer, intent(in) :: binary_id
            integer, intent(out) :: ierr
         end subroutine other_jdot_interface

         subroutine other_binary_wind_transfer_interface(binary_id, s_i, ierr)
            integer, intent(in) :: binary_id, s_i
            integer, intent(out) :: ierr
         end subroutine other_binary_wind_transfer_interface

         subroutine other_edot_interface(binary_id, ierr)
            integer, intent(in) :: binary_id
            integer, intent(out) :: ierr
         end subroutine other_edot_interface

         subroutine other_CE_init_interface(binary_id, restart, ierr)
            integer, intent(in) :: binary_id
            logical, intent(in) :: restart
            integer, intent(out) :: ierr
         end subroutine other_CE_init_interface

         subroutine other_CE_rlo_mdot_interface(binary_id, rlo_mdot, ierr)
            use const_def, only: dp
            integer, intent(in) :: binary_id
            real(dp), intent(out) :: rlo_mdot
            integer, intent(out) :: ierr
         end subroutine other_CE_rlo_mdot_interface

         integer function other_CE_binary_evolve_step_interface(binary_id)
            integer, intent(in) :: binary_id
         end function other_CE_binary_evolve_step_interface
         
         integer function other_CE_binary_finish_step_interface(binary_id)
            integer, intent(in) :: binary_id
         end function other_CE_binary_finish_step_interface
         
         integer function extras_binary_startup_interface(binary_id,restart,ierr)
            integer, intent(in) :: binary_id
            integer, intent(out) :: ierr
            logical,intent(in) :: restart    
         end function extras_binary_startup_interface
         
         integer function extras_binary_start_step_interface(binary_id, ierr)
            integer, intent(in) :: binary_id
            integer, intent(out) :: ierr
         end function extras_binary_start_step_interface
         
         integer function extras_binary_check_model_interface(binary_id)
            integer, intent(in) :: binary_id
         end function extras_binary_check_model_interface
         
         integer function extras_binary_finish_step_interface(binary_id)
            integer, intent(in) :: binary_id
         end function extras_binary_finish_step_interface

         subroutine extras_binary_after_evolve_interface(binary_id, ierr)
            integer, intent(in) :: binary_id
            integer, intent(out) :: ierr
         end subroutine extras_binary_after_evolve_interface

         subroutine other_binary_photo_write_interface(binary_id, iounit)
            integer, intent(in) :: binary_id, iounit
         end subroutine other_binary_photo_write_interface

         subroutine other_binary_photo_read_interface(binary_id, iounit, ierr)
            integer, intent(in) :: binary_id, iounit
            integer, intent(out) :: ierr
         end subroutine other_binary_photo_read_interface

         subroutine other_e2_interface(id, e2, ierr)
            use const_def, only: dp
            integer, intent(in) :: id
            real(dp),intent (out) :: e2
            integer, intent(out) :: ierr
         end subroutine other_e2_interface
               
      end interface

      type binary_job_controls
         include "binary_job_controls.inc"
      end type binary_job_controls

      type binary_info
         !binary id
         integer :: binary_id ! unique identifier for each binary_info instance
         logical :: in_use

         integer :: extra_binary_terminal_iounit

         type (binary_job_controls) :: job
         include 'binary_data.inc'
         include 'binary_controls.inc'
      end type

      logical :: have_initialized_binary_handles = .false.
      integer, parameter :: max_binary_handles = 10 ! this can be increased as necessary
      type (binary_info), target, save :: binary_handles(max_binary_handles)
         ! gfortran seems to require "save" here.  at least it did once upon a time.
      
      
      contains

      subroutine binary_ptr(binary_id, b, ierr)
         integer, intent(in) :: binary_id
         type (binary_info), pointer, intent(inout) :: b
         integer, intent(out) :: ierr
         call get_binary_ptr(binary_id, b, ierr)
      end subroutine binary_ptr


      subroutine get_binary_ptr(binary_id, b, ierr)
         integer, intent(in) :: binary_id
         type (binary_info), pointer :: b
         integer, intent(out) :: ierr         
         if (binary_id < 1 .or. binary_id > max_binary_handles) then
             ierr = -1
             return
         end if
         b => binary_handles(binary_id)
         ierr = 0
      end subroutine get_binary_ptr
      

      logical function is_donor(b, s)
         type (binary_info), pointer :: b
         type (star_info), pointer :: s
         is_donor = (s% id == b% d_i)
      end function is_donor
      
      subroutine init_binary_data
      
      end subroutine init_binary_data

      end module binary_def
