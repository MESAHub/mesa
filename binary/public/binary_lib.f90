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


module binary_lib
   
   use const_def
   
   implicit none

contains
   
   subroutine run1_binary(tst, &
      ! star extras
      extras_controls, &
      ! binary extras
      extras_binary_controls, &
      ierr, &
      inlist_fname_arg)
      
      use run_binary_support, only : do_run1_binary
      use binary_def, only : init_binary_data
      use star_def, only : star_info
      
      logical, intent(in) :: tst
      
      interface
         
         subroutine extras_controls(id, ierr)
            integer, intent(in) :: id
            integer, intent(out) :: ierr
         end subroutine extras_controls
         
         subroutine extras_binary_controls(binary_id, ierr)
            integer :: binary_id
            integer, intent(out) :: ierr
         end subroutine extras_binary_controls
      
      end interface
      
      integer, intent(out) :: ierr
      character (len = *) :: inlist_fname_arg
      optional inlist_fname_arg
      
      call init_binary_data
      
      call do_run1_binary(tst, &
         ! star extras
         extras_controls, &
         ! binary extras
         extras_binary_controls, &
         ierr, &
         inlist_fname_arg)
   
   end subroutine run1_binary
   
   subroutine binary_set_ignore_rlof_flag(binary_id, ignore_rlof_flag, ierr)
      use binary_utils, only : set_ignore_rlof_flag
      integer, intent(in) :: binary_id
      logical, intent(in) :: ignore_rlof_flag
      integer, intent(out) :: ierr
      
      ierr = 0
      call set_ignore_rlof_flag(binary_id, ignore_rlof_flag, ierr)
   end subroutine binary_set_ignore_rlof_flag
   
   subroutine binary_set_point_mass_i(binary_id, point_mass_i, ierr)
      use binary_utils, only : set_point_mass_i
      integer, intent(in) :: binary_id
      integer, intent(in) :: point_mass_i
      integer, intent(out) :: ierr
      
      ierr = 0
      call set_point_mass_i(binary_id, point_mass_i, ierr)
   end subroutine binary_set_point_mass_i
   
   subroutine binary_set_m1(binary_id, m1, ierr)
      use binary_utils, only : set_m1
      integer, intent(in) :: binary_id
      real(dp), intent(in) :: m1
      integer, intent(out) :: ierr
      
      ierr = 0
      call set_m1(binary_id, m1, ierr)
   end subroutine binary_set_m1
   
   subroutine binary_set_m2(binary_id, m2, ierr)
      use binary_utils, only : set_m2
      integer, intent(in) :: binary_id
      real(dp), intent(in) :: m2
      integer, intent(out) :: ierr
      
      ierr = 0
      call set_m2(binary_id, m2, ierr)
   end subroutine binary_set_m2
   
   subroutine binary_set_period_eccentricity(binary_id, period, eccentricity, ierr)
      use binary_utils, only : set_period_eccentricity
      integer, intent(in) :: binary_id
      real(dp), intent(in) :: period ! in seconds
      real(dp), intent(in) :: eccentricity
      integer, intent(out) :: ierr
      call set_period_eccentricity(binary_id, period, eccentricity, ierr)
   end subroutine binary_set_period_eccentricity
   
   subroutine binary_set_separation_eccentricity(binary_id, separation, eccentricity, ierr)
      use binary_utils, only : set_separation_eccentricity
      integer, intent(in) :: binary_id
      real(dp), intent(in) :: separation ! in cm
      real(dp), intent(in) :: eccentricity
      integer, intent(out) :: ierr
      call set_separation_eccentricity(binary_id, separation, eccentricity, ierr)
   
   end subroutine binary_set_separation_eccentricity
   
   real(dp) function binary_eval_rlobe(m1, m2, a)
      use binary_utils, only : eval_rlobe
      real(dp), intent(in) :: m1, m2, a
      ! Roche lobe size for star of mass m1 with a
      ! companion of mass m2 at separation a, according to
      ! the approximation of Eggleton 1983, apj 268:368-369
      binary_eval_rlobe = eval_rlobe(m1, m2, a)
   end function binary_eval_rlobe
   
   subroutine binary_eval_mdot_edd(binary_id, mdot_edd, mdot_edd_eta, ierr)
      use binary_mdot, only : eval_mdot_edd
      integer, intent(in) :: binary_id
      integer, intent(out) :: ierr
      real(dp), intent(out) :: mdot_edd, mdot_edd_eta
      call eval_mdot_edd(binary_id, mdot_edd, mdot_edd_eta, ierr)
   end subroutine binary_eval_mdot_edd
   
   subroutine binary_eval_accreted_material_j(binary_id, ierr)
      use binary_mdot, only : eval_accreted_material_j
      integer, intent(in) :: binary_id
      integer, intent(out) :: ierr
      call eval_accreted_material_j(binary_id, ierr)
   end subroutine binary_eval_accreted_material_j
   
   subroutine binary_eval_wind_xfer_fractions(binary_id, ierr)
      use binary_wind, only : eval_wind_xfer_fractions
      integer, intent(in) :: binary_id
      integer, intent(out) :: ierr
      call eval_wind_xfer_fractions(binary_id, ierr)
   end subroutine binary_eval_wind_xfer_fractions
   
   
   subroutine binary_get_control_namelist(binary_id, name, val, ierr)
      use binary_ctrls_io, only : get_binary_control
      use binary_def, only : binary_info, binary_ptr
      integer, intent(in) :: binary_id
      character(len = *), intent(in) :: name
      character(len = *), intent(out) :: val
      integer, intent(out) :: ierr
      type (binary_info), pointer :: b
      
      ierr = 0
      call binary_ptr(binary_id, b, ierr)
      if(ierr/=0) return
      call get_binary_control(b, name, val, ierr)
   
   end subroutine binary_get_control_namelist
   
   subroutine binary_set_control_namelist(binary_id, name, val, ierr)
      use binary_ctrls_io, only : set_binary_control
      use binary_def, only : binary_info, binary_ptr
      integer, intent(in) :: binary_id
      character(len = *), intent(in) :: name
      character(len = *), intent(in) :: val
      integer, intent(out) :: ierr
      type (binary_info), pointer :: b
      
      ierr = 0
      call binary_ptr(binary_id, b, ierr)
      if(ierr/=0) return
      call set_binary_control(b, name, val, ierr)
   
   end subroutine binary_set_control_namelist
   
   
   subroutine binary_get_star_job_namelist(binary_id, name, val, ierr)
      use binary_job_ctrls_io, only : get_binary_job
      use binary_def, only : binary_info, binary_ptr
      integer, intent(in) :: binary_id
      character(len = *), intent(in) :: name
      character(len = *), intent(out) :: val
      integer, intent(out) :: ierr
      type (binary_info), pointer :: b
      
      ierr = 0
      call binary_ptr(binary_id, b, ierr)
      if(ierr/=0) return
      call get_binary_job(b, name, val, ierr)
   
   end subroutine binary_get_star_job_namelist
   
   subroutine binary_set_star_job_namelist(binary_id, name, val, ierr)
      use binary_job_ctrls_io, only : set_binary_job
      use binary_def, only : binary_info, binary_ptr
      integer, intent(in) :: binary_id
      character(len = *), intent(in) :: name
      character(len = *), intent(in) :: val
      integer, intent(out) :: ierr
      type (binary_info), pointer :: b
      
      ierr = 0
      call binary_ptr(binary_id, b, ierr)
      if(ierr/=0) return
      call set_binary_job(b, name, val, ierr)
   
   end subroutine binary_set_star_job_namelist

end module binary_lib

