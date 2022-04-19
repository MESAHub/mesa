! ***********************************************************************
!
!   Copyright (C) 2011-2019  The MESA Team
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


      include 'overshoot_dbl_exp/overshoot_dbl_exp_def.inc'
      include 'timestep_limit/timestep_limit_def.inc'
      include 'other_winds/other_winds_def.inc'
      include 'xtrans_mesh_factor/xtrans_mesh_factor_def.inc'
      
      
      ! these routines are called by the standard run_star check_model
      contains

      include "test_suite_extras.inc"
      
      
      subroutine extras_controls(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         
         include 'overshoot_dbl_exp/overshoot_dbl_exp_extras_controls.inc'
         if (ierr /= 0) return
         include 'timestep_limit/timestep_limit_extras_controls.inc'
         if (ierr /= 0) return
         include 'other_winds/other_winds_extras_controls.inc'
         if (ierr /= 0) return
         include 'xtrans_mesh_factor/xtrans_mesh_factor_extras_controls.inc'
         if (ierr /= 0) return
         
         s% use_other_kap = .true.
         s% other_kap_get => my_kap_get
         
         s% eos_rq% use_other_eos_results = .true.
         s% eos_rq% other_eos_results => my_other_eos_results
         s% eos_rq% use_other_eos_component = .false.
         s% eos_rq% other_eos_frac => my_other_eos_frac
         s% eos_rq% other_eos_component => my_other_eos_component

         s% use_other_screening = .true.
         s% other_screening => my_screening

         s% use_other_rate_get = .true.
         s% other_rate_get => my_rate_get
         
         s% extras_startup => extras_startup
         s% extras_check_model => extras_check_model
         s% extras_finish_step => extras_finish_step
         s% extras_after_evolve => extras_after_evolve
         s% how_many_extra_history_columns => how_many_extra_history_columns
         s% data_for_extra_history_columns => data_for_extra_history_columns
         s% how_many_extra_profile_columns => how_many_extra_profile_columns
         s% data_for_extra_profile_columns => data_for_extra_profile_columns           
         
      end subroutine extras_controls


      include 'overshoot_dbl_exp/overshoot_dbl_exp.inc'
      include 'timestep_limit/timestep_limit.inc'
      include 'other_winds/other_winds.inc'
      include 'xtrans_mesh_factor/xtrans_mesh_factor.inc'
      
      
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
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
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
      

      ! eos hook routines

      subroutine my_other_eos_frac( &
              handle, &
              species, chem_id, net_iso, xa, &
              Rho, log10Rho, T, log10T, &
              frac, dfrac_dlogRho, dfrac_dlogT, ierr)

         ! INPUT
         use chem_def, only: num_chem_isos

         integer, intent(in) :: handle ! eos handle

         integer, intent(in) :: species
         integer, pointer :: chem_id(:) ! maps species to chem id
            ! index from 1 to species
            ! value is between 1 and num_chem_isos
         integer, pointer :: net_iso(:) ! maps chem id to species number
            ! index from 1 to num_chem_isos (defined in chem_def)
            ! value is 0 if the iso is not in the current net
            ! else is value between 1 and number of species in current net
         real(dp), intent(in) :: xa(:) ! mass fractions

         real(dp), intent(in) :: Rho, log10Rho ! the density
         real(dp), intent(in) :: T, log10T ! the temperature

         ! OUTPUT
         ! this routine must provide a fraction (in [0,1]) of the 'other' eos to use
         ! the remaining fraction (1-frac) will be provided by the standard MESA eos
         real(dp), intent(out) :: frac ! fraction of other_eos to use
         real(dp), intent(out) :: dfrac_dlogRho ! its partial derivative at constant T
         real(dp), intent(out) :: dfrac_dlogT   ! its partial derivative at constant Rho

         integer, intent(out) :: ierr ! 0 means AOK.

         ! this would use other_eos_component everywhere
         frac = 1d0
         dfrac_dlogRho = 0d0
         dfrac_dlogT = 0d0

      end subroutine my_other_eos_frac


      subroutine my_other_eos_component( &
              handle, &
              species, chem_id, net_iso, xa, &
              Rho, log10Rho, T, log10T, &
              res, d_dlnRho_c_T, d_dlnT_c_Rho, &
              d_dxa_c_TRho, ierr)

         use chem_def, only: num_chem_isos
         use eos_def
         use eos_lib

         ! INPUT

         integer, intent(in) :: handle ! eos handle
         integer, intent(in) :: species
         integer, pointer :: chem_id(:) ! maps species to chem id
         integer, pointer :: net_iso(:) ! maps chem id to species number
         real(dp), intent(in) :: xa(:) ! mass fractions
         real(dp), intent(in) :: Rho, log10Rho ! the density
         real(dp), intent(in) :: T, log10T ! the temperature
         real(dp), intent(inout) :: res(:) ! (num_eos_basic_results)
         real(dp), intent(inout) :: d_dlnRho_c_T(:) ! (num_eos_basic_results)
         real(dp), intent(inout) :: d_dlnT_c_Rho(:) ! (num_eos_basic_results)
         real(dp), intent(inout) :: d_dxa_c_TRho(:,:) ! (num_eos_d_dxa_results, species)

         integer, intent(out) :: ierr ! 0 means AOK.

         ! one must provide a complete set of EOS results
         ierr = -1

      end subroutine my_other_eos_component


      subroutine my_other_eos_results( &
              handle, &
              species, chem_id, net_iso, xa, &
              Rho, log10Rho, T, log10T, &
              res, d_dlnRho_c_T, d_dlnT_c_Rho, &
              d_dxa_c_TRho, ierr)

         use chem_def, only: num_chem_isos
         use eos_def
         use eos_lib

         ! INPUT

         integer, intent(in) :: handle ! eos handle
         integer, intent(in) :: species
         integer, pointer :: chem_id(:) ! maps species to chem id
         integer, pointer :: net_iso(:) ! maps chem id to species number
         real(dp), intent(in) :: xa(:) ! mass fractions
         real(dp), intent(in) :: Rho, log10Rho ! the density
         real(dp), intent(in) :: T, log10T ! the temperature
         real(dp), intent(inout) :: res(:) ! (num_eos_basic_results)
         real(dp), intent(inout) :: d_dlnRho_c_T(:) ! (num_eos_basic_results)
         real(dp), intent(inout) :: d_dlnT_c_Rho(:) ! (num_eos_basic_results)
         real(dp), intent(inout) :: d_dxa_c_TRho(:,:) ! (num_eos_d_dxa_results, species)

         integer, intent(out) :: ierr ! 0 means AOK.

         ! one can modify the existing eos results
         res(i_lnPgas) = res(i_lnPgas) + 0

      end subroutine my_other_eos_results
      
      
      subroutine my_kap_get( &
            id, k, handle, species, chem_id, net_iso, xa, &
            log10_rho, log10_T, &
            lnfree_e, d_lnfree_e_dlnRho, d_lnfree_e_dlnT, &
            eta, d_eta_dlnRho, d_eta_dlnT, &
            kap_fracs, kap, dln_kap_dlnRho, dln_kap_dlnT, dln_kap_dxa, ierr)

         use kap_def, only: num_kap_fracs
         use kap_lib
 
         ! INPUT
         integer, intent(in) :: id ! star id if available; 0 otherwise
         integer, intent(in) :: k ! cell number or 0 if not for a particular cell         
         integer, intent(in) :: handle ! from alloc_kap_handle
         integer, intent(in) :: species
         integer, pointer :: chem_id(:) ! maps species to chem id
            ! index from 1 to species
            ! value is between 1 and num_chem_isos         
         integer, pointer :: net_iso(:) ! maps chem id to species number
            ! index from 1 to num_chem_isos (defined in chem_def)
            ! value is 0 if the iso is not in the current net
            ! else is value between 1 and number of species in current net
         real(dp), intent(in) :: xa(:) ! mass fractions
         real(dp), intent(in) :: log10_rho ! density
         real(dp), intent(in) :: log10_T ! temperature
         real(dp), intent(in) :: lnfree_e, d_lnfree_e_dlnRho, d_lnfree_e_dlnT
            ! free_e := total combined number per nucleon of free electrons and positrons
         real(dp), intent(in) :: eta, d_eta_dlnRho, d_eta_dlnT
            ! eta := electron degeneracy parameter

         ! OUTPUT
         real(dp), intent(out) :: kap_fracs(num_kap_fracs)
         real(dp), intent(out) :: kap ! opacity
         real(dp), intent(out) :: dln_kap_dlnRho ! partial derivative at constant T
         real(dp), intent(out) :: dln_kap_dlnT   ! partial derivative at constant Rho
         real(dp), intent(out) :: dln_kap_dxa(:) ! partial derivative w.r.t. to species
         integer, intent(out) :: ierr ! 0 means AOK.
                  
         call kap_get( &
            handle, species, chem_id, net_iso, xa, log10_rho, log10_T, &
            lnfree_e, d_lnfree_e_dlnRho, d_lnfree_e_dlnT, &
            eta, d_eta_dlnRho, d_eta_dlnT, &
            kap_fracs, kap, dln_kap_dlnRho, dln_kap_dlnT, dln_kap_dxa, ierr)

      end subroutine my_kap_get


      subroutine my_screening(sc, z1, z2, a1, a2, screen, dscreendt, dscreendd, ierr)
         use rates_def
   
         implicit none
   
         type (Screen_Info) :: sc ! See rates_def
         ! This contains lots of useful things like temperature, density etc as well as some precomputed
         ! terms that are useful for screening calculations. The derived type is set in do_screen_set_context (screen.f90)
         real(dp),intent(in) ::    z1, z2      !< charge numbers of reactants
         real(dp),intent(in) ::    a1, a2     !< mass numbers of reactants
         real(dp),intent(out) ::   screen     !< on return, screening factor for this reaction
         real(dp),intent(out) ::   dscreendt     !< on return, temperature derivative of the screening factor
         real(dp),intent(out) ::   dscreendd    !< on return, density derivative of the screening factor
         integer, intent(out) ::   ierr
   
         screen = 1d0
         dscreendt = 0d0
         dscreendd = 0d0
         ierr = 0
      end subroutine my_screening


      subroutine my_rate_get(ir, temp, tf, raw_rate, ierr)
         use rates_def
         use rates_lib
         implicit none
   
         integer :: ir ! Rate id
         real(dp),intent(in) ::    temp      !< Temperature
         type (T_Factors) :: tf !< Various temperature factors
         real(dp),intent(inout) ::   raw_rate     !< Unscreened reaction_rate, note this will have the default mesa rate on entry
         integer, intent(out) ::   ierr
      
         ierr = 0
   
         if (.false. .and. trim(reaction_name(ir)) == 'r_he4_he4_he4_to_c12') then
            if(temp<1d8) then
               raw_rate = 0d0
            end if
         
         end if
   
      end subroutine my_rate_get


      end module run_star_extras
      
