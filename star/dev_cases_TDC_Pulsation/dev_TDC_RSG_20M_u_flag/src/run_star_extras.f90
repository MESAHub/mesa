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
      use chem_def
      use utils_lib
      use rates_def, only: i_rate
      use gyre_mesa_m


      use interp_1d_def, only: pm_work_size
      use interp_1d_lib, only: interp_pm, interp_values, interp_value

      implicit none
      
      include "test_suite_extras_def.inc"
      include 'run_star_extras_TDC_pulsation_defs.inc'

      logical :: dbg = .false.

      !!!!!!!!!!!!!!!!!!!!!!!!!
      ! These variables are loaded up from x_ctrl, x_integer_ctrl and x_logical_ctrl
      ! values specified on inlist_ppisn
      !!!!!!!!!!!!!!!!!!!!!!!!!


      logical :: in_inlist_pulses, remesh_for_envelope_model, turn_off_remesh, remove_core
      integer :: kick_model_number, timestep_drop_model_number, turn_off_remesh_model_number
      integer :: initial_model_number
      real(dp) :: max_dt_before_pulse, max_dt_during_pulse, core_T_for_cut

      contains

      include "test_suite_extras.inc"
      include 'run_star_extras_TDC_pulsation.inc'

      
      subroutine extras_controls(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         
         s% extras_startup => extras_startup
         s% extras_start_step => extras_start_step
         s% extras_check_model => extras_check_model
         s% extras_finish_step => extras_finish_step
         s% extras_after_evolve => extras_after_evolve
         s% how_many_extra_history_columns => how_many_extra_history_columns
         s% data_for_extra_history_columns => data_for_extra_history_columns
         s% how_many_extra_profile_columns => how_many_extra_profile_columns
         s% data_for_extra_profile_columns => data_for_extra_profile_columns

         ! pulsation info
         s% other_photo_write => photo_write
         s% other_photo_read => photo_read

         s% how_many_other_mesh_fcns => how_many_other_mesh_fcns
         s% other_mesh_fcn_data => RSP_mesh

         ! this is optional
         s% other_wind => brott_wind
         s% other_adjust_mdot => my_adjust_mdot
         s% other_before_struct_burn_mix => my_before_struct_burn_mix
         s% other_kap_get => my_other_kap_get
         s% other_alpha_mlt => alpha_mlt_routine

         ! store user provided options from the inlist

         in_inlist_pulses = s% x_logical_ctrl(22)
         max_dt_before_pulse = s% x_ctrl(17)
         max_dt_during_pulse = s% x_ctrl(18)
         remesh_for_envelope_model = s% x_logical_ctrl(23)
         turn_off_remesh = s% x_logical_ctrl(24)
         kick_model_number = s% x_ctrl(11)
         timestep_drop_model_number = s% x_ctrl(13)
         turn_off_remesh_model_number = s% x_ctrl(12)
         remove_core = s% x_logical_ctrl(25)
         core_T_for_cut = s% x_ctrl(14)
      end subroutine extras_controls


        subroutine alpha_mlt_routine(id, ierr)
        use chem_def, only: ih1
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
        do k=1,s% nz
        if (s% xa(h1,k) >= H_limit) then
        s% alpha_mlt(k) = alpha_H
        else
        s% alpha_mlt(k) = alpha_other
        end if
        !write(*,2) 'alpha_mlt', k, s% alpha_mlt(k),
        end do
        !stop
        end subroutine alpha_mlt_routine


      subroutine brott_wind(id, Lsurf, Msurf, Rsurf, Tsurf, X, Y, Z, w, ierr)
         use star_def
         integer, intent(in) :: id
         real(dp), intent(in) :: Lsurf, Msurf, Rsurf, Tsurf, X, Y, Z ! surface values (cgs)
         ! NOTE: surface is outermost cell. not necessarily at photosphere.
         ! NOTE: don't assume that vars are set at this point.
         ! so if you want values other than those given as args,
         ! you should use values from s% xh(:,:) and s% xa(:,:) only.
         ! rather than things like s% Teff or s% lnT(:) which have not been set yet.
         real(dp), intent(out) :: w ! wind in units of Msun/year (value is >= 0)
         integer, intent(out) :: ierr

         integer :: h1, he4
         real(dp) :: Xs, Ys, Z_div_Z_solar, Teff_jump, alfa, L1, M1, R1, T1, &
            vink_wind, nieu_wind, hamann_wind
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         L1 = Lsurf
         M1 = Msurf
         R1 = Rsurf
         T1 = Tsurf

         h1 = s% net_iso(ih1)
         he4 = s% net_iso(ihe4)
         Xs = s% xa(h1,1)
         Ys = s% xa(he4,1)
         ! Z=0.0142 is Z from Asplund et al. 2009
         Z_div_Z_solar = s% kap_rq% Zbase/0.0142d0
         ! use Vink et al 2001, eqns 14 and 15 to set "jump" temperature
         Teff_jump = 1d3*(61.2d0 + 2.59d0*(-13.636d0 + 0.889d0*log10(Z_div_Z_solar)))

         vink_wind = 0d0
         nieu_wind = 0d0
         hamann_wind = 0d0
         w = 0

         call eval_Vink_wind(vink_wind)
         call eval_Nieuwenhuijzen_wind(nieu_wind)
         call eval_Hamann_wind(hamann_wind)

         ! use 1/10 hamann
         hamann_wind = hamann_wind/10d0

         if (T1 < Teff_jump) then
            ! low T wind
            w = max(vink_wind, nieu_wind)
         else
            ! high T wind
            alfa = 0d0
            if (Xs > 0.7d0) then
               alfa = 1d0
            else if (Xs > 0.4d0 .and. Xs < 0.7d0) then
               alfa = (Xs - 0.4d0)/0.3d0
            end if
            w = alfa * vink_wind + (1d0-alfa) * hamann_wind
         end if

         ierr = 0

         contains

         subroutine eval_Vink_wind(w)
            real(dp), intent(inout) :: w
            real(dp) :: alfa, w1, w2, logMdot, dT, vinf_div_vesc

            ! alfa = 1 for hot side, = 0 for cool side
            if (T1 > 27500d0) then
               alfa = 1
            else if (T1 < 22500d0) then
               alfa = 0
            else
               dT = 100d0
               if (T1 > Teff_jump + dT) then
                  alfa = 1
               else if (T1 < Teff_jump - dT) then
                  alfa = 0
               else
                  alfa = (T1 - (Teff_jump - dT)) / (2*dT)
               end if
            end if
            
            if (alfa > 0) then ! eval hot side wind (eqn 24)
               vinf_div_vesc = 2.6d0 ! this is the hot side galactic value
               vinf_div_vesc = vinf_div_vesc*pow(Z_div_Z_solar,0.13d0) ! corrected for Z
               logMdot = &
                  - 6.697d0 &
                  + 2.194d0*log10(L1/Lsun/1d5) &
                  - 1.313d0*log10(M1/Msun/30) &
                  - 1.226d0*log10(vinf_div_vesc/2d0) &
                  + 0.933d0*log10(T1/4d4) &
                  - 10.92d0*pow2(log10(T1/4d4)) &
                  + 0.85d0*log10(Z_div_Z_solar)
               w1 = exp10(logMdot)
            else
               w1 = 0
            end if
            
            if (alfa < 1) then ! eval cool side wind (eqn 25)
               vinf_div_vesc = 1.3d0 ! this is the cool side galactic value
               vinf_div_vesc = vinf_div_vesc*pow(Z_div_Z_solar,0.13d0) ! corrected for Z
               logMdot = &
                  - 6.688d0 &
                  + 2.210d0*log10(L1/Lsun/1d5) &
                  - 1.339d0*log10(M1/Msun/30) &
                  - 1.601d0*log10(vinf_div_vesc/2d0) &
                  + 1.07d0*log10(T1/2d4) &
                  + 0.85d0*log10(Z_div_Z_solar)
               w2 = exp10(logMdot)
            else
               w2 = 0
            end if
            
            w = alfa*w1 + (1 - alfa)*w2
            
         end subroutine eval_Vink_wind

         subroutine eval_Nieuwenhuijzen_wind(w)
            ! Nieuwenhuijzen, H.; de Jager, C. 1990, A&A, 231, 134 (eqn 2)
            real(dp), intent(out) :: w
            real(dp) :: log10w
            include 'formats'
            log10w = -14.02d0 &
                     +1.24d0*log10(L1/Lsun) &
                     +0.16d0*log10(M1/Msun) &
                     +0.81d0*log10(R1/Rsun) &
                     +0.85d0*log10(Z_div_Z_solar)
            w = exp10(log10w)
         end subroutine eval_Nieuwenhuijzen_wind

         subroutine eval_Hamann_wind(w)
            ! Hamann, W.-R.; Koesterke, L.; Wessolowski, U. 1995, A&A, 299, 151
            real(dp), intent(out) :: w
            real(dp) :: log10w
            include 'formats'
            log10w = -11.95d0 &
                     +1.5d0*log10(L1/Lsun) &
                     -2.85d0*Xs &
                     + 0.85d0*log10(Z_div_Z_solar)
            w = exp10(log10w)
         end subroutine eval_Hamann_wind

      end subroutine brott_wind
      
      subroutine my_adjust_mdot(id, ierr)
         use star_def
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         real(dp) :: Lrad_div_Ledd
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         if (s% generations > 2) then
            write(*,*) "check mdots", s% mstar_dot, s% mstar_dot_old
            if (abs(s% mstar_dot) > 1.05d0*abs(s% mstar_dot_old)) then
               s% mstar_dot = 1.05d0*s% mstar_dot_old
            else if (abs(s% mstar_dot) < 0.95d0*abs(s% mstar_dot_old)) then
               s% mstar_dot = 0.95d0*s% mstar_dot_old
            end if
         end if
      end subroutine my_adjust_mdot


      subroutine my_other_kap_get( &
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

         type (star_info), pointer :: s
         real(dp) :: velocity
         real(dp) :: radius, logR
         real(dp) :: logT_alt, inv_diff
         real(dp) :: log_kap, alpha

         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
                  
         kap = 0; dln_kap_dlnRho = 0; dln_kap_dlnT = 0; dln_kap_dxa = 0
         velocity = 0
         radius = 0

         !if (k==1 .and. s% u_flag .and. .not. is_nan(s% lnR_start(1))) then !very surface cell can go mad, things are more stable if we fix opacity
         !   if (s% xh_start(s% i_u,1)>sqrt(2*s% cgrav(1)*s% m(1)/exp(s% lnR_start(1)))) then
         if (k==1 .and. s% u_flag) then !very surface cell can go mad, things are more stable if we fix opacity
            ! this is to support restarts, as xh_start and r_start are
            ! not properly set when model loads
            if (s% solver_iter > 0) then
               velocity = s% xh_start(s% i_u,1)
               radius = s% r_start(1)
            else
               velocity = s% xh(s% i_u,1)
               radius = s% r(1)
            end if
            if (velocity>sqrt(2*s% cgrav(1)*s% m(1)/radius)) then
               kap = 0.2d0*(1 + s% X(1))
               dln_kap_dlnRho = 0d0
               dln_kap_dlnT = 0d0
               return
            else
               call kap_get( &
                  s% kap_handle, species, chem_id, net_iso, xa, &
                  log10_rho, log10_T, lnfree_e, d_lnfree_e_dlnRho, d_lnfree_e_dlnT, &
                  eta, d_eta_dlnRho, d_eta_dlnT, &
                  kap_fracs, kap, dln_kap_dlnRho, dln_kap_dlnT, dln_kap_dxa, ierr)
            end if
         else
            call kap_get( &
               s% kap_handle, species, chem_id, net_iso, xa, &
               log10_rho, log10_T, lnfree_e, d_lnfree_e_dlnRho, d_lnfree_e_dlnT, &
               eta, d_eta_dlnRho, d_eta_dlnT, &
               kap_fracs, kap, dln_kap_dlnRho, dln_kap_dlnT, dln_kap_dxa, ierr)
         end if


      end subroutine my_other_kap_get

      subroutine extras_startup(id, restart, ierr)
         integer, intent(in) :: id
         logical, intent(in) :: restart
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         include 'formats'
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         call test_suite_startup(s, restart, ierr)
         call TDC_pulsation_extras_startup(id, restart, ierr)

         ! interestingly, if you do this instead of remove_center_by_temperature
         ! in the starjob section of the inlist, then the tau relaxation happens
         ! before the cut. Not sure which is better, but leaving like this for now.
         ! It turns out that this appears to be the better way and
         ! to do it, as this smooths the initialization of new atm BCs (if changed).
         if (.not. restart .and. in_inlist_pulses .and. remove_core) then
            call star_remove_center_by_temperature(id, core_T_for_cut, ierr)
          end if

         ! Initialize GYRE

         call init('gyre.in')

         ! Set constants

         call set_constant('G_GRAVITY', standard_cgrav)
         call set_constant('C_LIGHT', clight)
         call set_constant('A_RADIATION', crad)

         call set_constant('M_SUN', Msun)
         call set_constant('R_SUN', Rsun)
         call set_constant('L_SUN', Lsun)

         call set_constant('GYRE_DIR', TRIM(mesa_dir)//'/gyre/gyre')

         ! for rsp style mesh
         if (.not. restart .and. in_inlist_pulses ) then
         initial_model_number = s% model_number
         end if

         if (.not. restart .and. in_inlist_pulses .and. remesh_for_envelope_model) then
            initial_model_number = s% model_number

            call remesh_for_TDC_pulsation(id, ierr)
         end if
      end subroutine extras_startup
      
      
      subroutine extras_after_evolve(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         real(dp) :: dt
         character (len=strlen) :: test
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         call test_suite_after_evolve(s, ierr)

         if (.not. s% x_logical_ctrl(37)) return
         call final()
      end subroutine extras_after_evolve
      

      ! returns either keep_going, retry, or terminate.
      integer function extras_check_model(id)
         integer, intent(in) :: id
         integer :: ierr, k
         real(dp) :: max_v
         type (star_info), pointer :: s
         include 'formats'
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
         how_many_extra_history_columns = TDC_pulsation_how_many_extra_history_columns(id)
      end function how_many_extra_history_columns
      
      
      subroutine data_for_extra_history_columns(id, n, names, vals, ierr)
         integer, intent(in) :: id, n
         character (len=maxlen_history_column_name) :: names(n)
         real(dp) :: vals(n), v_esc
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer :: k, k0
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         call TDC_pulsation_data_for_extra_history_columns(id, n, names, vals, ierr)
      end subroutine data_for_extra_history_columns

      
      integer function how_many_extra_profile_columns(id)
         use star_def, only: star_info
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
!         how_many_extra_profile_columns = TDC_pulsation_how_many_extra_profile_columns(id)
         how_many_extra_profile_columns = 1!how_many_extra_profile_columns +1
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
!         call TDC_pulsation_data_for_extra_profile_columns(id, n, nz, names, vals, ierr)
         names(1) = 'Frad_div_cUrad'
         do k = 1, s% nz
         vals(k,1) = ((s% L(k) - s% L_conv(k)) / (4._dp*pi*pow2(s%r(k)))) &
            /(clight * s% Prad(k) *3._dp)
         end do
      end subroutine data_for_extra_profile_columns
      

      integer function extras_start_step(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         include 'formats'
         extras_start_step = terminate
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         !this is used to ensure we read the right inlist options
         s% use_other_before_struct_burn_mix = .true.

         ! we want to ignore T gradient equation for a few steps after remesh
         if (s% model_number < initial_model_number + 100 .and. in_inlist_pulses) then
            s% convergence_ignore_equL_residuals = .true.
         else if (in_inlist_pulses) then
            s% convergence_ignore_equL_residuals = .false.
         end if

         if (s% model_number == kick_model_number .and. in_inlist_pulses &
            .and. s% x_logical_ctrl(5))then

            ! if v= 0, turn on v so we can kick
            if (.not. s% v_flag .or. .not. s% u_flag) then
               call star_set_v_flag(id, .true., ierr)
            end if

            call gyre_in_mesa_extras_set_velocities(s,.false.,ierr)
            write(*,*) 'kick'
            write(*,*) 'kick'
            write(*,*) 'kick'
            write(*,*) 'kick'
            write(*,*) 'kick'

         end if

         call my_before_struct_burn_mix(s% id, s% dt, extras_start_step)

        ! add stopping condition for testing.
        if ((.not. in_inlist_pulses) .and. s% center_h1 < 1d-1) then
            s% Teff_lower_limit = exp10(3.5563d0) ! joyce et al. 2020.
            !write(*,*) 'stopping because have reached ~ 3600 K Teff for Betelgeuse'
        else
            s% Teff_lower_limit = -1d99
        end if

         extras_start_step = keep_going
      end function extras_start_step
   
      subroutine my_before_struct_burn_mix(id, dt, res)
         use const_def, only: dp
         use star_def
         integer, intent(in) :: id
         real(dp), intent(in) :: dt
         integer, intent(out) :: res ! keep_going, redo, retry, terminate
         real(dp) :: power_photo, v_esc
         integer :: ierr, k
         type (star_info), pointer :: s
         include 'formats'
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

!            !s% use_atm_PT_at_center_of_surface_cell = .false.
!            s% use_momentum_outer_BC = .false. ! offset_P_to_center_cell = .true.
!            s% use_compression_outer_BC = .false. ! offset_P_to_center_cell = .true.
!            s% use_zero_Pgas_outer_BC = .true.
!            s% atm_option = 'T_tau'
!            s% atm_T_tau_relation = 'Eddington'
!            s% atm_T_tau_opacity = 'fixed'
!            s% tau_factor = 1d-3
!            s% Pextra_factor = 1d0
!            s% force_tau_factor = 1d-3
!            s% delta_lgL_limit = 0.25d0
!            !s% delta_lgTeff_limit = 1d-2!0.25d0
!            s% delta_lgL_limit_L_min = 1d99!-100
!            s% delta_lgL_limit_L_min = 1d99!-100
!
!            !s% atm_T_tau_errtol = 1d-12
!            !s% atm_T_tau_max_iters = 500
  
         if (in_inlist_pulses) then
            if (s% model_number > timestep_drop_model_number )then
                 s% max_timestep = max_dt_during_pulse
            else
                 s% max_timestep = max_dt_before_pulse
            end if

         ! time step control on pulsations
         if (period > 0d0 .and. period/s% max_timestep < 600 .and. &
         s% model_number > timestep_drop_model_number) then
             s% max_timestep = period/600d0
         end if

            if (s% model_number > turn_off_remesh_model_number .and. turn_off_remesh )then
               s% okay_to_remesh = .false.
         !               if ((s% model_number == turn_off_remesh_model_number + 1) &
         !                     .and. .not. remesh_for_envelope_model )then
         !                  do k =1,s%nz
         !                     if (s%lnT(k) >= log(2d6)) then
         !                        exit
         !                     end if
         !                  end do
         !                  s% mesh_min_k_old_for_split = k
         !               end if
         !               write (*,*) 's% mesh_min_k_old_for_split', s% mesh_min_k_old_for_split
            end if
         end if
         
         ! adjust convection near surface to prevent crazy issues:

!         do k = 1, s% nz
!            if (s% tau(k) > 2d0/3d0) then
!               exit
!            end if
!         end do
!
!         s% max_q_for_convection_with_hydro_on = s% q(k-1)

         ! reading inlists can turn this flag off for some reason
         s% use_other_before_struct_burn_mix = .true.

         res = keep_going
      end subroutine my_before_struct_burn_mix
      
      subroutine null_binary_controls(id, binary_id, ierr)
         integer, intent(in) :: id, binary_id
         integer, intent(out) :: ierr
         ierr = 0
      end subroutine null_binary_controls

      ! returns either keep_going or terminate.
      integer function extras_finish_step(id)
         use run_star_support
         use math_lib
         integer, intent(in) :: id
         integer :: ierr,k
         real(dp) :: max_vel_inside, vesc_for_cell, vesc_surf !check_avg_v_div_vesc
         type (star_info), pointer :: s
         include 'formats'
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         extras_finish_step = keep_going
         extras_finish_step = TDC_pulsation_extras_finish_step(id)

!         if (.not. s% x_logical_ctrl(37)) return
!         extras_finish_step = gyre_in_mesa_extras_finish_step(id)

         if (extras_finish_step == terminate) s% termination_code = t_extras_finish_step

      end function extras_finish_step




      ! here is an example that adds a mesh function for log(opacity)
      subroutine how_many_other_mesh_fcns(id, n)
         integer, intent(in) :: id
         integer, intent(out) :: n
         n = 1
      end subroutine how_many_other_mesh_fcns


!      subroutine RSP_mesh( &
!         id, nfcns, names, gval_is_xa_function, vals1, ierr)
!         use star_def
!         use math_lib
!         use const_def
!         integer, intent(in) :: id
!         integer, intent(in) :: nfcns
!         character(len=*) :: names(:)
!         logical, intent(out) :: gval_is_xa_function(:)  ! (nfcns)
!         real(dp), pointer :: vals1(:)  ! =(nz, nfcns)
!         integer, intent(out) :: ierr
!         integer :: nz, k
!         real(dp), pointer :: vals(:, :)
!         real(dp) :: weight1 = 1d6!1d4
!         real(dp) :: weight2 = 8d5!1d4
!         real(dp) :: weight3 = 0d0
!         real(dp) :: logT_anchor1, logT_anchor2, logT_anchor3, lmid, delta, ell
!         integer :: k_anchor1, k_anchor2
!
!         type(star_info), pointer :: s
!         ierr = 0
!         call star_ptr(id, s, ierr)
!         if (ierr /= 0) return
!         names(1) = 'RSP_function'
!         gval_is_xa_function(1) = .false.
!         nz = s%nz
!         vals(1:nz, 1:nfcns) => vals1(1:nz*nfcns)
!
!         logT_anchor1 = log(11d3)!log(11d3)
!         logT_anchor2 = log(20d3)!log(11d3)
!         logT_anchor3 = log(30d3)
!
!         lmid  = 0.5d0*(logT_anchor2+logT_anchor3)
!         delta = (logT_anchor3 - logT_anchor2)
!
!         k_anchor1 = 0
!         k_anchor2 = 0
!
!      !         do k = 1, nz
!      !            if (s% lnT(k) < logT_anchor1) then
!      !               vals(k, 1) = weight1*(s% m(1) - s% m(k))/Msun
!      !               k_anchor = k
!      !               !write (*,*) "k", k ,"dm", vals(k, 1)
!      !            else if (s% lnT(k) < logT_anchor2 .and. s% lnT(k) >= logT_anchor1) then
!      !               vals(k, 1) = weight2*(s% m(k_anchor) - s% m(k))/Msun
!      !               !write (*,*) "k", k ,"dm", vals(k, 1)
!      !            else
!      !               vals(k, 1) = weight3*(s% m(1)/Msun - s% m(k)/Msun)
!      !            end if
!      !         end do
!
!
!         do k = 1, nz
!            ell = s%lnT(k)
!            if (s% lnT(k) <= logT_anchor1) then
!               vals(k,1) = weight1*(s%m(1) - s%m(k))/Msun
!               k_anchor1 = k
!            else if (s% lnT(k) <= logT_anchor2 .and. s% lnT(k) >= logT_anchor1) then
!               vals(k,1) = weight2*(s%m(1) - s%m(k))/Msun
!               k_anchor2 = k
!            else if (s% lnT(k) < logT_anchor3) then
!               ! smooth taper doqn to 0.
!!               vals(k,1) = vals(k-1,1)/2d0
!               vals(k,1) = (0.5d0*weight2*(1d0 - tanh( (ell - lmid)/delta )) &
!                  ) * ( (s%m(k_anchor2) - s%m(k)) / Msun )
!            end if
!         end do
!
!      end subroutine RSP_mesh

subroutine RSP_mesh( &
   id, nfcns, names, gval_is_xa_function, vals1, ierr)
   use star_def
   use math_lib
   use const_def
   integer, intent(in) :: id
   integer, intent(in) :: nfcns
   character(len=*) :: names(:)
   logical, intent(out) :: gval_is_xa_function(:)  ! (nfcns)
   real(dp), pointer :: vals1(:)  ! =(nz, nfcns)
   integer, intent(out) :: ierr
   integer :: nz, k
   real(dp), pointer :: vals(:, :)
   real(dp) :: weight1 = 1d5!1d4
   real(dp) :: weight2 = 1d5!1d4
   real(dp) :: weight3 = 0d0
   real(dp) :: logT_anchor1, logT_anchor2, logT_anchor3, lmid, delta, ell
   integer :: k_anchor1, k_anchor2

   type(star_info), pointer :: s
   ierr = 0
   call star_ptr(id, s, ierr)
   if (ierr /= 0) return
   names(1) = 'RSP_function'
   gval_is_xa_function(1) = .false.
   nz = s%nz
   vals(1:nz, 1:nfcns) => vals1(1:nz*nfcns)

   logT_anchor1 = log(11d3)!log(11d3)
   logT_anchor2 = log(20d3)!log(11d3)
   logT_anchor3 = log(30d3)

   lmid  = 0.5d0*(logT_anchor2+logT_anchor3)
   delta = (logT_anchor3 - logT_anchor2)/10d0

   k_anchor1 = 0
   k_anchor2 = 0

!         do k = 1, nz
!            if (s% lnT(k) < logT_anchor1) then
!               vals(k, 1) = weight1*(s% m(1) - s% m(k))/Msun
!               k_anchor = k
!               !write (*,*) "k", k ,"dm", vals(k, 1)
!            else if (s% lnT(k) < logT_anchor2 .and. s% lnT(k) >= logT_anchor1) then
!               vals(k, 1) = weight2*(s% m(k_anchor) - s% m(k))/Msun
!               !write (*,*) "k", k ,"dm", vals(k, 1)
!            else
!               vals(k, 1) = weight3*(s% m(1)/Msun - s% m(k)/Msun)
!            end if
!         end do


   do k = 1, nz
      ell = s%lnT(k)
      if (s% lnT(k) < logT_anchor1) then
         ! core weighting
         vals(k,1) = weight1*(s%m(1) - s%m(k))/Msun
         k_anchor1 = k
      else if (s% lnT(k) < logT_anchor2 .and. s% lnT(k) >= logT_anchor1) then
         ! envelope weighting
         vals(k,1) = weight2*(s%m(k_anchor1) - s%m(k))/Msun
         k_anchor2 = k
      else
         ! smooth taper from weight2 → 0
         vals(k,1) = ( &
            0.5d0*weight2*(1d0 - tanh( (ell - lmid)/delta )) &
            ) * ( (s%m(k_anchor2) - s%m(k)) / Msun )
      end if
   end do

end subroutine RSP_mesh


      subroutine photo_write(id, iounit)
         integer, intent(in) :: id, iounit
         call TDC_pulsation_photo_write(id, iounit)
      end subroutine photo_write


      subroutine photo_read(id, iounit, ierr)
         integer, intent(in) :: id, iounit
         integer, intent(out) :: ierr
         call TDC_pulsation_photo_read(id, iounit, ierr)
      end subroutine photo_read

      end module run_star_extras
      
