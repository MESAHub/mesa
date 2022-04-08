! ***********************************************************************
!
!   Copyright (C) 2010-2019  Bill Paxton & The MESA Team
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
      use chem_def !! maybe taking up uncessary space but whatever
      

      implicit none
      real(dp) :: spotf, spotx, PB_i, Teff_local, sb_sigma

      include "test_suite_extras_def.inc"

      contains

      include "test_suite_extras.inc"
   

      ! these routines are called by the standard run_star check_model
     ! contains
      
      subroutine extras_controls(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         
         ! this is the place to set any procedure pointers you want to change
         ! e.g., other_wind, other_mixing, other_energy  (see star_data.inc)

         ! the extras functions in this file will not be called
         ! unless you set their function pointers as done below.
         ! otherwise we use a null_ version which does nothing (except warn).


         spotf = s% x_ctrl(1)
         spotx = s% x_ctrl(2)


         s% extras_startup => extras_startup
         s% extras_start_step => extras_start_step
         s% extras_check_model => extras_check_model
         s% extras_finish_step => extras_finish_step
         s% extras_after_evolve => extras_after_evolve
         s% how_many_extra_history_columns => how_many_extra_history_columns
         s% data_for_extra_history_columns => data_for_extra_history_columns
         s% how_many_extra_profile_columns => how_many_extra_profile_columns
         s% data_for_extra_profile_columns => data_for_extra_profile_columns  

         s% how_many_extra_history_header_items => how_many_extra_history_header_items
         s% data_for_extra_history_header_items => data_for_extra_history_header_items
         s% how_many_extra_profile_header_items => how_many_extra_profile_header_items
         s% data_for_extra_profile_header_items => data_for_extra_profile_header_items


         if (s% use_other_surface_PT) then
            s% other_surface_PT => starspot_tweak_PT
            !s% other_surface_PT => starspot_tweak_PT_Joyce
         end if


         if (s% use_other_mlt_results) then
           s% other_mlt_results => YREC_spots_other_mlt_results
         end if
        

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
      

      integer function extras_start_step(id)
         integer, intent(in) :: id
         integer :: ierr
         real(dp) :: mu_ideal_gas, R2, R_gas_constant !! pi can be deleted
         type (star_info), pointer :: s
         real(dp) :: power_he_burn, power_c_burn, power_neutrinos, &
         center_h1, center_he4, ocz_top_mass, ocz_bot_mass, &
         ocz_top_radius, ocz_bot_radius!, mass_difference_prev !! no!!
         integer :: nz, j, i, k, k_ocz_top, k_ocz_bot, n_conv_bdy
!         include 'formats'
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_start_step = 0

         !! to use constants as defined in MESA, import const_def and call
         !! variables by their names there

         ! set PB_i here and then other_mlt will know about it, because 
         ! quantities set here are carried over the course of a timestep, 
         ! AKA all newton iterations
         
         sb_sigma = boltz_sigma !5.67d-5 !!cgs units
         mu_ideal_gas = s% mu(1)  !1.00794d0 ! for hydrogen, 1 gram per mole
         !write(*,*) 'MESA def, my def: ', s% mu(1), mu_ideal_gas  
         R_gas_constant = cgas !8.314d7 ! cgs units

         R2 = pow2(s%R(1))
         Teff_local = pow( s%L(1)/(4.0*pi*sb_sigma*R2), 0.25d0)
         PB_i = (R_gas_constant* s%rho(1)/mu_ideal_gas) * (1.0 - spotx) * Teff_local

      end function extras_start_step


!-----------------------------------------------------------------------------------------
!
! YREC routines
!
!----------------------------------------------------------------------------------------
      subroutine YREC_spots_other_mlt_results(id, k, MLT_option, &  ! NOTE: k=0 is a valid arg
            r, L, T, P, opacity, rho, chiRho, chiT, Cp, gradr, grada, scale_height, &
            iso, XH1, cgrav, m, gradL_composition_term, mixing_length_alpha, &
            alpha_semiconvection, thermohaline_coeff, &
            mixing_type, gradT, Y_face, conv_vel, D, Gamma, ierr)
         use const_def, only: dp
         use auto_diff
         integer, intent(in) :: id
         integer, intent(in) :: k
         ! integer, intent(out) :: ierr
         ! type (star_info), pointer :: s
         ! integer :: nz, j
         character (len=*), intent(in) :: MLT_option
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
        type(auto_diff_real_star_order1) :: spotx_of_r !, spotx4
        type(auto_diff_real_star_order1) :: gradr_spot
            !ierr = 0
        ! ------------------------------ 10/26/21
        type (star_info), pointer :: s
        ierr = 0
        call star_ptr(id, s, ierr)
        if (ierr /= 0) return
         
        !------------------------------
         !if (s% star_age >= 10d0) then
         if (.not. s% doing_relax) then
            spotx_of_r = (P - PB_i)/P 
            gradr_spot = gradr/( spotf*pow( spotx_of_r, 4d0) + 1d0 - spotf)
         else
            gradr_spot = gradr
         end if

         call star_mlt_results(id, k, MLT_option, &
            r, L, T, P, opacity, rho, chiRho, chiT, Cp, gradr_spot, grada, scale_height, &
            iso, XH1, cgrav, m, gradL_composition_term, mixing_length_alpha, &
            alpha_semiconvection, thermohaline_coeff, &
            mixing_type, gradT, Y_face, conv_vel, D, Gamma, ierr)
      end subroutine YREC_spots_other_mlt_results


      subroutine starspot_tweak_PT(id, skip_partials, &
            lnT_surf, dlnT_dL, dlnT_dlnR, dlnT_dlnM, dlnT_dlnkap, &
            lnP_surf, dlnP_dL, dlnP_dlnR, dlnP_dlnM, dlnP_dlnkap, ierr)

            use const_def, only: dp

         
            integer, intent(in) :: id
            logical, intent(in) :: skip_partials

            logical :: need_atm_Psurf, need_atm_Tsurf
         
            real(dp), intent(out) :: &
                  lnT_surf, dlnT_dL, dlnT_dlnR, dlnT_dlnM, dlnT_dlnkap, &
                  lnP_surf, dlnP_dL, dlnP_dlnR, dlnP_dlnM, dlnP_dlnkap
            integer, intent(out) :: ierr

            ! For my tweaks
            real(dp) ::  alp !, spotf, spotx

            ! Call the stock get_surf_PT
            type (star_info), pointer :: s
            real(dp) :: L_init

            include 'formats'

            ierr = 0
            call star_ptr(id, s, ierr)
            if (ierr /= 0) return

            need_atm_Psurf = .true.
            need_atm_Tsurf = .true. 
            sb_sigma = boltz_sigma
            
            alp = 1d0 - spotf + spotf*spotx*spotx*spotx*spotx

            ! This is the surface-average value for luminosity
            L_init = s% L(1)

            ! Set the surface L to the unspotted, ambient L
            s% L(1) = s% L(1) / alp

            ! Now, set the Teff. Used in atm table lookup to set boundary conditions
            s% Teff = pow(s% L(1)/(4._dp*pi*s% r(1)*s% r(1)*sb_sigma), 0.25_dp)

            ! Set everything with Lamb.
            call star_get_surf_PT(id, skip_partials, need_atm_Psurf, need_atm_Tsurf, &
                  lnT_surf, dlnT_dL, dlnT_dlnR, dlnT_dlnM, dlnT_dlnkap, &
                  lnP_surf, dlnP_dL, dlnP_dlnR, dlnP_dlnM, dlnP_dlnkap, &
                  ierr)

            s% Teff = pow(L_init/(4._dp*pi*s% r(1)*s% r(1)*sb_sigma), 0.25_dp)
            s% L(1) = L_init

      end subroutine starspot_tweak_PT

!------------------------------------------------------------------------------
!
! end YREC routines 
!
!------------------------------------------------------------------------------

      ! returns either keep_going, retry, or terminate.
      integer function extras_check_model(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_check_model = keep_going     


         if (.false. .and. s% star_mass_h1 < 0.35d0) then
            ! stop when star hydrogen mass drops to specified level
            extras_check_model = terminate
            write(*, *) 'have reached desired hydrogen mass'
            return
         end if

         ! by default, indicate where (in the code) MESA terminated
         if (extras_check_model == terminate) s% termination_code = t_extras_check_model
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
         
         ! note: do NOT add the extras names to history_columns.list
         ! the history_columns.list is only for the built-in history column options.
         ! it must not include the new column names you are adding here.
         

      end subroutine data_for_extra_history_columns

      
      integer function how_many_extra_profile_columns(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_profile_columns = 0
      end function how_many_extra_profile_columns
      
      
      subroutine data_for_extra_profile_columns(id, n, nz, names, vals, ierr)
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


      integer function how_many_extra_history_header_items(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_history_header_items = 0
      end function how_many_extra_history_header_items


      subroutine data_for_extra_history_header_items(id, n, names, vals, ierr)
         integer, intent(in) :: id, n
         character (len=maxlen_history_column_name) :: names(n)
         real(dp) :: vals(n)
         type(star_info), pointer :: s
         integer, intent(out) :: ierr
         ierr = 0
         call star_ptr(id,s,ierr)
         if(ierr/=0) return
      end subroutine data_for_extra_history_header_items


      integer function how_many_extra_profile_header_items(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_profile_header_items = 0
      end function how_many_extra_profile_header_items


      subroutine data_for_extra_profile_header_items(id, n, names, vals, ierr)
         integer, intent(in) :: id, n
         character (len=maxlen_profile_column_name) :: names(n)
         real(dp) :: vals(n)
         type(star_info), pointer :: s
         integer, intent(out) :: ierr
         ierr = 0
         call star_ptr(id,s,ierr)
         if(ierr/=0) return
      end subroutine data_for_extra_profile_header_items


      ! returns either keep_going or terminate.
      ! note: cannot request retry; extras_check_model can do that.
      integer function extras_finish_step(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_finish_step = keep_going

         ! to save a profile, 
            ! s% need_to_save_profiles_now = .true.
         ! to update the star log,
            ! s% need_to_update_history_now = .true.

         ! see extras_check_model for information about custom termination codes
         ! by default, indicate where (in the code) MESA terminated
         if (extras_finish_step == terminate) s% termination_code = t_extras_finish_step
      end function extras_finish_step
      
      
      subroutine extras_after_evolve(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         call test_suite_after_evolve(s, ierr)
      end subroutine extras_after_evolve

      ! subroutine extras_after_evolve(id, ierr)
      !    integer, intent(in) :: id
      !    integer, intent(out) :: ierr
      !    type (star_info), pointer :: s
      !    real(dp) :: dt, eEmg
      !    integer :: k, k0
         
      !    ierr = 0
      !    call star_ptr(id, s, ierr)
      !    if (ierr /= 0) return
         
         

      !    ! Find a zone at mass coordinate 0.5
      !    do k=1,s%nz
      !       if(s% m(k)/msun < 0.5d0) then
      !          k0 = k
      !          exit
      !       end if
      !    end do
      !    eEmg = qe * s% E_field(k0)/(amu * s% g_field_element_diffusion(k0))
         
      !    write(*,*) 'Core eE/mg = ', eEmg
      !    if (eEmg > 1.95d0 .and. eEmg < 2.05d0) then
      !       write(*,*) 'passed test for electric field in the core'
      !    else
      !       write(*,*) 'failed test for electric field in the core'
      !    end if         
         
      !    call test_suite_after_evolve(s, ierr)
      ! end subroutine extras_after_evolve



      end module run_star_extras
      
