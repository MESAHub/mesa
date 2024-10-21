! ***********************************************************************
!
!   Copyright (C) 2010-2024  The MESA Team
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
   use turb_def

   implicit none

   type(th_info_t), allocatable, save :: th_info(:)

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

      s% other_mlt_results => th_other_mlt_results
      
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
      integer :: k, k_cntr, k_surf

      logical :: okay
      type (star_info), pointer :: s
      character (len=strlen) :: test
      include 'formats'
      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return
      okay = .true.

      call test_suite_after_evolve(s, ierr)
      write(*,'(A)')
      call check('log total_angular_momentum', safe_log10(s% total_angular_momentum), 45d0, 55d0)
      call check('log center_omega', safe_log10(s% center_omega), -4d0, -2d0)
      call check('log he_core_omega', safe_log10(s% he_core_omega), -6d0, -3d0)
      call check('surface j_rot', safe_log10(s% j_rot(1)),  5d0, 25d0)
      call check('surface v_rot', s% omega(1)*s% r(1)*1d-5, 0d0, 1d0)

      k_cntr = 0
      k_surf = 0
      do k = s% nz, 1, -1
         if (s% m(k) > 0.24d0*Msun .and. k_cntr == 0) k_cntr = k
         if (s% m(k) > 0.25d0*Msun .and. k_surf == 0) k_surf = k
      end do

      if (k_surf >=1 .and. k_cntr <= s% nz) then
         write(*,'(A)')
         write(*,1) 'avg near 0.245 Msun'
         call check('logT', avg_val(s% lnT)/ln10, 7d0, 7.5d0)
         call check('logRho', avg_val(s% lnd)/ln10, 0.5d0, 2.0d0)
         call check('log j_rot', safe_log10(avg_val(s% j_rot)), 10d0, 20d0)
         call check('D_ST', safe_log10(avg_val(s% D_ST)), 1d0, 8d0)
         call check('nu_ST', safe_log10(avg_val(s% nu_ST)), 4.0d0, 9.0d0)
         !call check('dynamo_B_r', safe_log10(avg_val(s% dynamo_B_r)), 0d0, 2d0)
         !call check('dynamo_B_phi', safe_log10(avg_val(s% dynamo_B_phi)), 3d0, 7d0)
         write(*,'(A)')
         if (okay) write(*,'(a)') 'all values are within tolerances'
      end if
      write(*,'(A)')


   contains

      real(dp) function avg_val(v)
         real(dp) :: v(:)
         avg_val = dot_product(v(k_surf:k_cntr), s% dq(k_surf:k_cntr)) / sum(s% dq(k_surf:k_cntr))
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
      include 'formats'

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
      how_many_extra_profile_columns = 20
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

      ! Extract data from th_info

      names = [ &
         'K_therm  ', & ! Thermal conductivity
         'K_T      ', & ! Thermal diffusivity
         'K_C      ', & ! Chemical diffusivity
         'nu       ', & ! Viscosity
         'Pr       ', & ! Prandtl number
         'tau      ', & ! Chemical diffusivity ratio
         'R_0      ', & ! Density ratio
         'r        ', & ! Reduced density ratio
         'H_B      ', & ! Lorentz force coefficient
         'Pm       ', & ! Magnetic Prandtl number
         'D_B      ', & ! Magnetic diffusivity
         'lam_hat  ', & ! Growth rate of fastest-growing fingering mode
         'l2_hat   ', & ! Horizontal wavenumber squared of fastest-growing fingering mode
         'sigma_max', & ! Growth rate of fastest-growing parasitic mode
         'k_z_max  ', & ! Vertical wavenumber of fastest-growing parasitic mode
         'w        ', & ! Saturation flow speed
         'w_HG19   ', & ! Saturation flow speed in HG19 treatment
         'w_FRG24  ', & ! Saturation flow speed in FRG24 treatment
         'Nu_C     ', & ! Compositional Nusselt number
         'D_thrm   ' &  ! Effective thermohaline mixing diffusivity 
         ]

      do k = 1, s%nz

         vals(k,:) = [ &
            th_info(k)%K_therm, &
            th_info(k)%K_T, &
            th_info(k)%K_C, &
            th_info(k)%nu, &
            th_info(k)%Pr, &
            th_info(k)%tau, &
            th_info(k)%R_0, &
            th_info(k)%r, &
            th_info(k)%H_B, &
            th_info(k)%Pm, &
            th_info(k)%D_B, &
            th_info(k)%lam_hat, &
            th_info(k)%l2_hat, &
            th_info(k)%sigma_max, &
            th_info(k)%k_z_max, &
            th_info(k)%w, &
            th_info(k)%w_HG19, &
            th_info(k)%w_FRG24, &
            th_info(k)%Nu_C, &
            th_info(k)%D_thrm &
            ]

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


   ! adds hook in to store th_info data
   subroutine th_other_mlt_results(id, k, MLT_option, &  ! NOTE: k=0 is a valid arg
      r, L, T, P, opacity, rho, chiRho, chiT, Cp, gradr, grada, scale_height, &
      iso, XH1, cgrav, m, gradL_composition_term, mixing_length_alpha, &
      alpha_semiconvection, thermohaline_coeff, &
      mixing_type, gradT, Y_face, conv_vel, D, Gamma, ierr)
      use const_def, only: dp
      use auto_diff
      use star_def
      use turb_def
      integer, intent(in) :: id
      integer, intent(in) :: k
      character(len=*), intent(in) :: MLT_option
      type(auto_diff_real_star_order1), intent(in) :: &
         r, L, T, P, opacity, rho, chiRho, chiT, Cp, gradr, grada, scale_height
      integer, intent(in) :: iso
      real(dp), intent(in) :: &
         XH1, cgrav, m, gradL_composition_term, &
         mixing_length_alpha, alpha_semiconvection, thermohaline_coeff
      integer, intent(out) :: mixing_type
      type(auto_diff_real_star_order1), intent(out) :: &
         gradT, Y_face, conv_vel, D, Gamma
      integer, intent(out) :: ierr

      type(star_info), pointer :: s
      
      ierr = 0

      ! Get the star pointer
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return

      ! If necessary allocate/resize the th_info array (throwing away all existing data!)
      !$OMP CRITICAL
      if (ALLOCATED(th_info)) then
         if (SIZE(th_info) /= s% nz) then
            deallocate(th_info)
            allocate(th_info(s% nz))
         end if
      else
         allocate(th_info(s% nz))
      end if
      !$OMP END CRITICAL

      ! Dispatch to star_mlt_results
      call star_mlt_results(id, k, MLT_option, &  ! NOTE: k=0 is a valid arg
         r, L, T, P, opacity, rho, chiRho, chiT, Cp, gradr, grada, scale_height, &
         iso, XH1, cgrav, m, gradL_composition_term, mixing_length_alpha, &
         alpha_semiconvection, thermohaline_coeff, &
         mixing_type, gradT, Y_face, conv_vel, D, Gamma, ierr, th_info(k))
         
   end subroutine th_other_mlt_results

end module run_star_extras

