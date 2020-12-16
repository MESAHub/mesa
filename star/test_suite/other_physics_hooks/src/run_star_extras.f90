! ***********************************************************************
!
!   Copyright (C) 2011-2019  Bill Paxton & The MESA Team
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


      include 'convective_bdy/convective_bdy_def.inc'
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
         
         include 'convective_bdy/convective_bdy_extras_controls.inc'
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
         
         s% use_other_eos = .true.
         s% other_eosDT_get => my_eosDT_get
         s% other_eosDT_get_T => my_eosDT_get_T
         s% other_eosDT_get_Rho => my_eosDT_get_Rho
         s% other_eosPT_get => my_eosPT_get
         s% other_eosPT_get_T => my_eosPT_get_T
         s% other_eosPT_get_Pgas => my_eosPT_get_Pgas
         s% other_eosPT_get_Pgas_for_Rho => my_eosPT_get_Pgas_for_Rho
         
         s% extras_startup => extras_startup
         s% extras_check_model => extras_check_model
         s% extras_finish_step => extras_finish_step
         s% extras_after_evolve => extras_after_evolve
         s% how_many_extra_history_columns => how_many_extra_history_columns
         s% data_for_extra_history_columns => data_for_extra_history_columns
         s% how_many_extra_profile_columns => how_many_extra_profile_columns
         s% data_for_extra_profile_columns => data_for_extra_profile_columns           
         
      end subroutine extras_controls


      include 'convective_bdy/convective_bdy.inc'
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
      
      





      
      subroutine my_eosDT_get( &
              id, k, handle, Z, X, abar, zbar, & 
              species, chem_id, net_iso, xa, &
              Rho, log10Rho, T, log10T, & 
              res, d_dlnRho_c_T, d_dlnT_c_Rho, &
              d_dabar_c_TRho, d_dzbar_c_TRho, ierr)

         use chem_def, only: num_chem_isos
         use eos_def
         use eos_lib
         
         ! INPUT
         
         integer, intent(in) :: id ! star id if available; 0 otherwise
         integer, intent(in) :: k ! cell number or 0 if not for a particular cell         
         integer, intent(in) :: handle ! eos handle
         real(dp), intent(in) :: Z ! the metals mass fraction
         real(dp), intent(in) :: X ! the hydrogen mass fraction
         real(dp), intent(in) :: abar, zbar
         integer, intent(in) :: species
         integer, pointer :: chem_id(:) ! maps species to chem id
         integer, pointer :: net_iso(:) ! maps chem id to species number
         real(dp), intent(in) :: xa(:) ! mass fractions
         real(dp), intent(in) :: Rho, log10Rho ! the density
         real(dp), intent(in) :: T, log10T ! the temperature
         real(dp), intent(inout) :: res(:) ! (num_eos_basic_results)
         real(dp), intent(inout) :: d_dlnRho_c_T(:) ! (num_eos_basic_results) 
         real(dp), intent(inout) :: d_dlnT_c_Rho(:) ! (num_eos_basic_results) 
         real(dp), intent(inout) :: d_dabar_c_TRho(:) ! (num_eos_basic_results) 
         real(dp), intent(inout) :: d_dzbar_c_TRho(:) ! (num_eos_basic_results) 
         
         integer, intent(out) :: ierr ! 0 means AOK.
         
         call eosDT_get( &
              handle, Z, X, abar, zbar, & 
              species, chem_id, net_iso, xa, &
              Rho, log10Rho, T, log10T, & 
              res, d_dlnRho_c_T, d_dlnT_c_Rho, &
              d_dabar_c_TRho, d_dzbar_c_TRho, ierr)
         
      end subroutine my_eosDT_get
      
      
      ! the following routine uses gas pressure and temperature as input variables
      subroutine my_eosPT_get(&
               id, k, handle, Z, X, abar, zbar, &
               species, chem_id, net_iso, xa,&
               Pgas, log10Pgas, T, log10T, &
               Rho, log10Rho, dlnRho_dlnPgas_const_T, dlnRho_dlnT_const_Pgas, &
               res, d_dlnRho_const_T, d_dlnT_const_Rho, &
               d_dabar_c_TRho, d_dzbar_c_TRho, ierr)

         use eos_def
         use eos_lib
         use chem_def, only: num_chem_isos

         ! INPUT
         
         integer, intent(in) :: id ! star id if available; 0 otherwise
         integer, intent(in) :: k ! cell number or 0 if not for a particular cell         
         integer, intent(in) :: handle

         real(dp), intent(in) :: Z ! the metals mass fraction
         real(dp), intent(in) :: X ! the hydrogen mass fraction
            
         real(dp), intent(in) :: abar, zbar
         
         integer, intent(in) :: species
         integer, pointer :: chem_id(:) ! maps species to chem id
            ! index from 1 to species
            ! value is between 1 and num_chem_isos         
         integer, pointer :: net_iso(:) ! maps chem id to species number
            ! index from 1 to num_chem_isos (defined in chem_def)
            ! value is 0 if the iso is not in the current net
            ! else is value between 1 and number of species in current net
         real(dp), intent(in) :: xa(:) ! mass fractions
         
         real(dp), intent(in) :: Pgas, log10Pgas ! the gas pressure
            ! provide both if you have them.  else pass one and set the other to arg_not_provided
            ! "arg_not_provided" is defined in mesa const_def
            
         real(dp), intent(in) :: T, log10T ! the temperature
            ! provide both if you have them.  else pass one and set the other to arg_not_provided
                     
         ! OUTPUT
         
         real(dp), intent(out) :: Rho, log10Rho ! density
         real(dp), intent(out) :: dlnRho_dlnPgas_const_T
         real(dp), intent(out) :: dlnRho_dlnT_const_Pgas
         real(dp), intent(inout) :: res(:) ! (num_eos_basic_results)
         ! partial derivatives of the basic results wrt lnd and lnT
         real(dp), intent(inout) :: d_dlnRho_const_T(:) ! (num_eos_basic_results) 
         ! d_dlnRho_const_T(i) = d(res(i))/dlnd|T
         real(dp), intent(inout) :: d_dlnT_const_Rho(:) ! (num_eos_basic_results) 
         ! d_dlnT_const_Rho(i) = d(res(i))/dlnT|Rho
         real(dp), intent(inout) :: d_dabar_c_TRho(:) ! (num_eos_basic_results) 
         real(dp), intent(inout) :: d_dzbar_c_TRho(:) ! (num_eos_basic_results) 
         
         integer, intent(out) :: ierr ! 0 means AOK.
         
         call eosPT_get( &
            handle, Z, X, abar, zbar, &
            species, chem_id, net_iso, xa,&
            Pgas, log10Pgas, T, log10T, &
            Rho, log10Rho, dlnRho_dlnPgas_const_T, dlnRho_dlnT_const_Pgas, &
            res, d_dlnRho_const_T, d_dlnT_const_Rho, &
            d_dabar_c_TRho, d_dzbar_c_TRho, ierr)
         
      end subroutine my_eosPT_get
      
      

      ! eosDT search routines.
      
      subroutine my_eosDT_get_T( &
               id, k, handle, Z, X, abar, zbar, &
               species, chem_id, net_iso, xa, &
               logRho, which_other, other_value, &
               logT_tol, other_tol, max_iter, logT_guess, & 
               logT_bnd1, logT_bnd2, other_at_bnd1, other_at_bnd2, &
               logT_result, res, d_dlnRho_c_T, d_dlnT_c_Rho, &
               d_dabar_c_TRho, d_dzbar_c_TRho, &
               eos_calls, ierr)
     
         ! finds log10 T given values for density and 'other', and initial guess for temperature.
         ! does up to max_iter attempts until logT changes by less than tol.
         
         ! 'other' can be any of the basic result variables for the eos
         ! specify 'which_other' by means of the definitions in eos_def (e.g., i_lnE)
         
         use eos_def
         use eos_lib

         integer, intent(in) :: id ! star id if available; 0 otherwise
         integer, intent(in) :: k ! cell number or 0 if not for a particular cell         
         integer, intent(in) :: handle

         real(dp), intent(in) :: Z ! the metals mass fraction
         real(dp), intent(in) :: X ! the hydrogen mass fraction
            
         real(dp), intent(in) :: abar, zbar
         
         integer, intent(in) :: species
         integer, pointer :: chem_id(:) ! maps species to chem id
            ! index from 1 to species
            ! value is between 1 and num_chem_isos         
         integer, pointer :: net_iso(:) ! maps chem id to species number
            ! index from 1 to num_chem_isos (defined in chem_def)
            ! value is 0 if the iso is not in the current net
            ! else is value between 1 and number of species in current net
         real(dp), intent(in) :: xa(:) ! mass fractions
         
         real(dp), intent(in) :: logRho ! log10 of density
         integer, intent(in) :: which_other ! from eos_def.  e.g., i_lnE
         real(dp), intent(in) :: other_value ! desired value for the other variable
         real(dp), intent(in) :: other_tol
         
         real(dp), intent(in) :: logT_tol
         integer, intent(in) :: max_iter ! max number of iterations        

         real(dp), intent(in) :: logT_guess ! log10 of temperature
         real(dp), intent(in) :: logT_bnd1, logT_bnd2 ! bounds for logT
            ! if don't know bounds, just set to arg_not_provided (defined in const_def)
         real(dp), intent(in) :: other_at_bnd1, other_at_bnd2 ! values at bounds
            ! if don't know these values, just set to arg_not_provided (defined in const_def)

         real(dp), intent(out) :: logT_result ! log10 of temperature
         real(dp), intent(inout) :: res(:) ! (num_eos_basic_results)
         real(dp), intent(inout) :: d_dlnRho_c_T(:) ! (num_eos_basic_results)
         real(dp), intent(inout) :: d_dlnT_c_Rho(:) ! (num_eos_basic_results)
         real(dp), intent(inout) :: d_dabar_c_TRho(:) ! (num_eos_basic_results) 
         real(dp), intent(inout) :: d_dzbar_c_TRho(:) ! (num_eos_basic_results) 
         
         integer, intent(out) :: eos_calls
         integer, intent(out) :: ierr ! 0 means AOK.
         
         call eosDT_get_T( &
            handle, Z, X, abar, zbar, &
            species, chem_id, net_iso, xa, &
            logRho, which_other, other_value, &
            logT_tol, other_tol, max_iter, logT_guess, & 
            logT_bnd1, logT_bnd2, other_at_bnd1, other_at_bnd2, &
            logT_result, res, d_dlnRho_c_T, d_dlnT_c_Rho, &
            d_dabar_c_TRho, d_dzbar_c_TRho, &
            eos_calls, ierr)
         

      end subroutine my_eosDT_get_T
      

      subroutine my_eosDT_get_Rho( &
               id, k, handle, Z, X, abar, zbar, &
               species, chem_id, net_iso, xa, &
               logT, which_other, other_value, &
               logRho_tol, other_tol, max_iter, logRho_guess,  &
               logRho_bnd1, logRho_bnd2, other_at_bnd1, other_at_bnd2, &
               logRho_result, res, d_dlnRho_c_T, d_dlnT_c_Rho, &
               d_dabar_c_TRho, d_dzbar_c_TRho, eos_calls, ierr)
     
         ! finds log10 Rho given values for temperature and 'other', and initial guess for density.
         ! does up to max_iter attempts until logRho changes by less than tol.
         
         ! 'other' can be any of the basic result variables for the eos
         ! specify 'which_other' by means of the definitions in eos_def (e.g., i_lnE)
         
         use chem_def, only: num_chem_isos
         use eos_def
         use eos_lib
         
         integer, intent(in) :: id ! star id if available; 0 otherwise
         integer, intent(in) :: k ! cell number or 0 if not for a particular cell         
         integer, intent(in) :: handle

         real(dp), intent(in) :: Z ! the metals mass fraction
         real(dp), intent(in) :: X ! the hydrogen mass fraction
            
         real(dp), intent(in) :: abar, zbar
         
         integer, intent(in) :: species
         integer, pointer :: chem_id(:) ! maps species to chem id
            ! index from 1 to species
            ! value is between 1 and num_chem_isos         
         integer, pointer :: net_iso(:) ! maps chem id to species number
            ! index from 1 to num_chem_isos (defined in chem_def)
            ! value is 0 if the iso is not in the current net
            ! else is value between 1 and number of species in current net
         real(dp), intent(in) :: xa(:) ! mass fractions
         
         real(dp), intent(in) :: logT ! log10 of temperature

         integer, intent(in) :: which_other ! from eos_def.  e.g., i_lnE
         real(dp), intent(in) :: other_value ! desired value for the other variable
         real(dp), intent(in) :: other_tol
         
         real(dp), intent(in) :: logRho_tol

         integer, intent(in) :: max_iter ! max number of Newton iterations        

         real(dp), intent(in) :: logRho_guess ! log10 of density
         real(dp), intent(in) :: logRho_bnd1, logRho_bnd2 ! bounds for logRho
            ! if don't know bounds, just set to arg_not_provided (defined in const_def)
         real(dp), intent(in) :: other_at_bnd1, other_at_bnd2 ! values at bounds
            ! if don't know these values, just set to arg_not_provided (defined in const_def)

         real(dp), intent(out) :: logRho_result ! log10 of density

         real(dp), intent(inout) :: res(:) ! (num_eos_basic_results)
         real(dp), intent(inout) :: d_dlnRho_c_T(:) ! (num_eos_basic_results)
         real(dp), intent(inout) :: d_dlnT_c_Rho(:) ! (num_eos_basic_results)
         real(dp), intent(inout) :: d_dabar_c_TRho(:) ! (num_eos_basic_results) 
         real(dp), intent(inout) :: d_dzbar_c_TRho(:) ! (num_eos_basic_results) 

         integer, intent(out) :: eos_calls
         integer, intent(out) :: ierr ! 0 means AOK.
         
         call eosDT_get_Rho( &
            handle, Z, X, abar, zbar, &
            species, chem_id, net_iso, xa, &
            logT, which_other, other_value, &
            logRho_tol, other_tol, max_iter, logRho_guess,  &
            logRho_bnd1, logRho_bnd2, other_at_bnd1, other_at_bnd2, &
            logRho_result, res, d_dlnRho_c_T, d_dlnT_c_Rho, &
            d_dabar_c_TRho, d_dzbar_c_TRho, eos_calls, ierr)
         

      end subroutine my_eosDT_get_Rho
      
      
      
      ! eosPT search routines
      
      subroutine my_eosPT_get_T( &
               id, k, handle, Z, X, abar, zbar, &
               species, chem_id, net_iso, xa,&
               logPgas, which_other, other_value,&
               logT_tol, other_tol, max_iter, logT_guess, &
               logT_bnd1, logT_bnd2, other_at_bnd1, other_at_bnd2,&
               logT_result, Rho, log10Rho, dlnRho_dlnPgas_const_T, dlnRho_dlnT_const_Pgas, &
               res, d_dlnRho_const_T, d_dlnT_const_Rho, &
               d_dabar_c_TRho, d_dzbar_c_TRho, eos_calls, ierr)
     
         ! finds log10 T given values for gas pressure and 'other',
         ! and initial guess for temperature.
         ! does up to max_iter attempts until logT changes by less than tol.
         
         ! 'other' can be any of the basic result variables for the eos
         ! specify 'which_other' by means of the definitions in eos_def (e.g., i_lnE)
         
         use chem_def, only: num_chem_isos
         use eos_def
         use eos_lib
         
         integer, intent(in) :: id ! star id if available; 0 otherwise
         integer, intent(in) :: k ! cell number or 0 if not for a particular cell         
         integer, intent(in) :: handle

         real(dp), intent(in) :: Z ! the metals mass fraction
         real(dp), intent(in) :: X ! the hydrogen mass fraction
            
         real(dp), intent(in) :: abar, zbar
         
         integer, intent(in) :: species
         integer, pointer :: chem_id(:) ! maps species to chem id
            ! index from 1 to species
            ! value is between 1 and num_chem_isos         
         integer, pointer :: net_iso(:) ! maps chem id to species number
            ! index from 1 to num_chem_isos (defined in chem_def)
            ! value is 0 if the iso is not in the current net
            ! else is value between 1 and number of species in current net
         real(dp), intent(in) :: xa(:) ! mass fractions
         
         real(dp), intent(in) :: logPgas ! log10 of gas pressure
         integer, intent(in) :: which_other ! from eos_def.  e.g., i_lnE
         real(dp), intent(in) :: other_value ! desired value for the other variable
         real(dp), intent(in) :: other_tol
         
         real(dp), intent(in) :: logT_tol
         integer, intent(in) :: max_iter ! max number of iterations        

         real(dp), intent(in) :: logT_guess ! log10 of temperature
         real(dp), intent(in) :: logT_bnd1, logT_bnd2 ! bounds for logT
            ! if don't know bounds, just set to arg_not_provided (defined in const_def)
         real(dp), intent(in) :: other_at_bnd1, other_at_bnd2 ! values at bounds
            ! if don't know these values, just set to arg_not_provided (defined in const_def)

         real(dp), intent(out) :: logT_result ! log10 of temperature
         real(dp), intent(out) :: Rho, log10Rho ! density
         real(dp), intent(out) :: dlnRho_dlnPgas_const_T
         real(dp), intent(out) :: dlnRho_dlnT_const_Pgas

         real(dp), intent(inout) :: res(:) ! (num_eos_basic_results)
         real(dp), intent(inout) :: d_dlnRho_const_T(:) ! (num_eos_basic_results)
         real(dp), intent(inout) :: d_dlnT_const_Rho(:) ! (num_eos_basic_results)
         real(dp), intent(inout) :: d_dabar_c_TRho(:) ! (num_eos_basic_results) 
         real(dp), intent(inout) :: d_dzbar_c_TRho(:) ! (num_eos_basic_results) 
         
         integer, intent(out) :: eos_calls
         integer, intent(out) :: ierr ! 0 means AOK.
         
         call eosPT_get_T( &
            handle, Z, X, abar, zbar, &
            species, chem_id, net_iso, xa,&
            logPgas, which_other, other_value,&
            logT_tol, other_tol, max_iter, logT_guess, &
            logT_bnd1, logT_bnd2, other_at_bnd1, other_at_bnd2,&
            logT_result, Rho, log10Rho, dlnRho_dlnPgas_const_T, dlnRho_dlnT_const_Pgas, &
            res, d_dlnRho_const_T, d_dlnT_const_Rho, &
            d_dabar_c_TRho, d_dzbar_c_TRho, &
            eos_calls, ierr)
         
               
      end subroutine my_eosPT_get_T
      

      subroutine my_eosPT_get_Pgas(&
               id, k, handle, Z, X, abar, zbar, &
               species, chem_id, net_iso, xa,&
               logT, which_other, other_value,&
               logPgas_tol, other_tol, max_iter, logPgas_guess, &
               logPgas_bnd1, logPgas_bnd2, other_at_bnd1, other_at_bnd2,&
               logPgas_result, Rho, log10Rho, dlnRho_dlnPgas_const_T, dlnRho_dlnT_const_Pgas, &
               res, d_dlnRho_const_T, d_dlnT_const_Rho, &
               d_dabar_c_TRho, d_dzbar_c_TRho, eos_calls, ierr)
     
         ! finds log10 Pgas given values for temperature and 'other', and initial guess for gas pressure.
         ! does up to max_iter attempts until logPgas changes by less than tol.
         
         ! 'other' can be any of the basic result variables for the eos
         ! specify 'which_other' by means of the definitions in eos_def (e.g., i_lnE)
         
         use chem_def, only: num_chem_isos
         use eos_def
         use eos_lib
         
         integer, intent(in) :: id ! star id if available; 0 otherwise
         integer, intent(in) :: k ! cell number or 0 if not for a particular cell         
         integer, intent(in) :: handle

         real(dp), intent(in) :: Z ! the metals mass fraction
         real(dp), intent(in) :: X ! the hydrogen mass fraction
            
         real(dp), intent(in) :: abar, zbar
         
         integer, intent(in) :: species
         integer, pointer :: chem_id(:) ! maps species to chem id
            ! index from 1 to species
            ! value is between 1 and num_chem_isos         
         integer, pointer :: net_iso(:) ! maps chem id to species number
            ! index from 1 to num_chem_isos (defined in chem_def)
            ! value is 0 if the iso is not in the current net
            ! else is value between 1 and number of species in current net
         real(dp), intent(in) :: xa(:) ! mass fractions
         
         real(dp), intent(in) :: logT ! log10 of temperature

         integer, intent(in) :: which_other ! from eos_def.  e.g., i_lnE
         real(dp), intent(in) :: other_value ! desired value for the other variable
         real(dp), intent(in) :: other_tol
         
         real(dp), intent(in) :: logPgas_tol

         integer, intent(in) :: max_iter ! max number of Newton iterations        

         real(dp), intent(in) :: logPgas_guess ! log10 of gas pressure
         real(dp), intent(in) :: logPgas_bnd1, logPgas_bnd2 ! bounds for logPgas
            ! if don't know bounds, just set to arg_not_provided (defined in const_def)
         real(dp), intent(in) :: other_at_bnd1, other_at_bnd2 ! values at bounds
            ! if don't know these values, just set to arg_not_provided (defined in const_def)

         real(dp), intent(out) :: logPgas_result ! log10 of gas pressure
         real(dp), intent(out) :: Rho, log10Rho ! density
         real(dp), intent(out) :: dlnRho_dlnPgas_const_T
         real(dp), intent(out) :: dlnRho_dlnT_const_Pgas

         real(dp), intent(inout) :: res(:) ! (num_eos_basic_results)
         real(dp), intent(inout) :: d_dlnRho_const_T(:) ! (num_eos_basic_results)
         real(dp), intent(inout) :: d_dlnT_const_Rho(:) ! (num_eos_basic_results)
         real(dp), intent(inout) :: d_dabar_c_TRho(:) ! (num_eos_basic_results) 
         real(dp), intent(inout) :: d_dzbar_c_TRho(:) ! (num_eos_basic_results) 

         integer, intent(out) :: eos_calls
         integer, intent(out) :: ierr ! 0 means AOK.
         
         call eosPT_get_Pgas(&
            handle, Z, X, abar, zbar, &
            species, chem_id, net_iso, xa,&
            logT, which_other, other_value,&
            logPgas_tol, other_tol, max_iter, logPgas_guess, &
            logPgas_bnd1, logPgas_bnd2, other_at_bnd1, other_at_bnd2,&
            logPgas_result, Rho, log10Rho, dlnRho_dlnPgas_const_T, dlnRho_dlnT_const_Pgas, &
            res, d_dlnRho_const_T, d_dlnT_const_Rho, &
            d_dabar_c_TRho, d_dzbar_c_TRho, eos_calls, ierr)
         
         
      end subroutine my_eosPT_get_Pgas
      

      subroutine my_eosPT_get_Pgas_for_Rho(&
               id, k, handle, Z, X, abar, zbar, &
               species, chem_id, net_iso, xa,&
               logT, logRho_want,&
               logPgas_tol, logRho_tol, max_iter, logPgas_guess, &
               logPgas_bnd1, logPgas_bnd2, logRho_at_bnd1, logRho_at_bnd2,&
               logPgas_result, Rho, logRho, dlnRho_dlnPgas_const_T, dlnRho_dlnT_const_Pgas, &
               res, d_dlnRho_const_T, d_dlnT_const_Rho, &
               d_dabar_c_TRho, d_dzbar_c_TRho, eos_calls, ierr)
     
         ! finds log10 Pgas given values for temperature and density, and initial guess for gas pressure.
         ! does up to max_iter attempts until logPgas changes by less than tol.
         
         use chem_def, only: num_chem_isos         
         use eos_def
         use eos_lib
         
         integer, intent(in) :: id ! star id if available; 0 otherwise
         integer, intent(in) :: k ! cell number or 0 if not for a particular cell         
         integer, intent(in) :: handle

         real(dp), intent(in) :: Z ! the metals mass fraction
         real(dp), intent(in) :: X ! the hydrogen mass fraction
            
         real(dp), intent(in) :: abar, zbar
         
         integer, intent(in) :: species
         integer, pointer :: chem_id(:) ! maps species to chem id
            ! index from 1 to species
            ! value is between 1 and num_chem_isos         
         integer, pointer :: net_iso(:) ! maps chem id to species number
            ! index from 1 to num_chem_isos (defined in chem_def)
            ! value is 0 if the iso is not in the current net
            ! else is value between 1 and number of species in current net
         real(dp), intent(in) :: xa(:) ! mass fractions
         
         real(dp), intent(in) :: logT ! log10 of temperature

         real(dp), intent(in) :: logRho_want ! log10 of desired density
         real(dp), intent(in) :: logRho_tol
         
         real(dp), intent(in) :: logPgas_tol

         integer, intent(in) :: max_iter ! max number of Newton iterations        

         real(dp), intent(in) :: logPgas_guess ! log10 of gas pressure
         real(dp), intent(in) :: logPgas_bnd1, logPgas_bnd2 ! bounds for logPgas
            ! if don't know bounds, just set to arg_not_provided (defined in const_def)
         real(dp), intent(in) :: logRho_at_bnd1, logRho_at_bnd2 ! values at bounds
            ! if don't know these values, just set to arg_not_provided (defined in const_def)

         real(dp), intent(out) :: logPgas_result ! log10 of gas pressure
         real(dp), intent(out) :: Rho, logRho ! density corresponding to logPgas_result
         real(dp), intent(out) :: dlnRho_dlnPgas_const_T
         real(dp), intent(out) :: dlnRho_dlnT_const_Pgas

         real(dp), intent(inout) :: res(:) ! (num_eos_basic_results)
         real(dp), intent(inout) :: d_dlnRho_const_T(:) ! (num_eos_basic_results)
         real(dp), intent(inout) :: d_dlnT_const_Rho(:) ! (num_eos_basic_results)
         real(dp), intent(inout) :: d_dabar_c_TRho(:) ! (num_eos_basic_results) 
         real(dp), intent(inout) :: d_dzbar_c_TRho(:) ! (num_eos_basic_results) 

         integer, intent(out) :: eos_calls
         integer, intent(out) :: ierr ! 0 means AOK.

         call eosPT_get_Pgas_for_Rho( &
            handle, Z, X, abar, zbar, &
            species, chem_id, net_iso, xa,&
            logT, logRho_want,&
            logPgas_tol, logRho_tol, max_iter, logPgas_guess, &
            logPgas_bnd1, logPgas_bnd2, logRho_at_bnd1, logRho_at_bnd2,&
            logPgas_result, Rho, logRho, dlnRho_dlnPgas_const_T, dlnRho_dlnT_const_Pgas, &
            res, d_dlnRho_const_T, d_dlnT_const_Rho, &
            d_dabar_c_TRho, d_dzbar_c_TRho, eos_calls, ierr)
         
         
      end subroutine my_eosPT_get_Pgas_for_Rho


      subroutine my_kap_get( &
            id, k, handle, zbar, X, Z, XC, XN, XO, XNe, &
            log10_rho, log10_T, species, chem_id, net_iso, xa, &
            lnfree_e, d_lnfree_e_dlnRho, d_lnfree_e_dlnT, &
            frac_Type2, kap, dln_kap_dlnRho, dln_kap_dlnT, ierr)
         
         use kap_lib
 
         ! INPUT
         integer, intent(in) :: id ! star id if available; 0 otherwise
         integer, intent(in) :: k ! cell number or 0 if not for a particular cell         
         integer, intent(in) :: handle ! from alloc_kap_handle
         real(dp), intent(in) :: zbar ! average ion charge
         real(dp), intent(in) :: X, Z, XC, XN, XO, XNe ! composition    
         real(dp), intent(in) :: log10_rho ! density
         real(dp), intent(in) :: log10_T ! temperature
         real(dp), intent(in) :: lnfree_e, d_lnfree_e_dlnRho, d_lnfree_e_dlnT
            ! free_e := total combined number per nucleon of free electrons and positrons

         integer, intent(in) :: species
         integer, pointer :: chem_id(:) ! maps species to chem id
            ! index from 1 to species
            ! value is between 1 and num_chem_isos         
         integer, pointer :: net_iso(:) ! maps chem id to species number
            ! index from 1 to num_chem_isos (defined in chem_def)
            ! value is 0 if the iso is not in the current net
            ! else is value between 1 and number of species in current net
         real(dp), intent(in) :: xa(:) ! mass fractions
         
         ! OUTPUT
         real(dp), intent(out) :: frac_Type2
         real(dp), intent(out) :: kap ! opacity
         real(dp), intent(out) :: dln_kap_dlnRho ! partial derivative at constant T
         real(dp), intent(out) :: dln_kap_dlnT   ! partial derivative at constant Rho
         integer, intent(out) :: ierr ! 0 means AOK.
                  
         call kap_get( &
            handle, zbar, X, Z, XC, XN, XO, XNe, log10_rho, log10_T, &
            lnfree_e, d_lnfree_e_dlnRho, d_lnfree_e_dlnT, &
            frac_Type2, kap, dln_kap_dlnRho, dln_kap_dlnT, ierr)

      end subroutine my_kap_get


      end module run_star_extras
      
