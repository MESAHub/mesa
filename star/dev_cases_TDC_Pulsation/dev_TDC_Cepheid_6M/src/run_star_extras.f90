! ***********************************************************************
!
!   Copyright (C) 2010-2025  The MESA Team
!
!   This program is free software: you can redistribute it and/or modify
!   it under the terms of the GNU Lesser General Public License
!   as published by the Free Software Foundation,
!   either version 3 of the License, or (at your option) any later version.
!
!   This program is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!   See the GNU Lesser General Public License for more details.
!
!   You should have received a copy of the GNU Lesser General Public License
!   along with this program. If not, see <https://www.gnu.org/licenses/>.
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
   ! values specified on inlist_common, inlist_pulses
      !!!!!!!!!!!!!!!!!!!!!!!!!

   logical :: in_inlist_pulses, remesh_for_envelope_model, turn_off_remesh
   integer :: kick_model_number, timestep_drop_model_number, turn_off_remesh_model_number
   integer :: initial_model_number
   real(dp) :: max_dt_before_pulse, max_dt_during_pulse

contains

   include "test_suite_extras.inc"
   include 'run_star_extras_TDC_pulsation.inc'

   subroutine extras_controls(id, ierr)
      integer, intent(in) :: id
      integer, intent(out) :: ierr
      type(star_info), pointer :: s
      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return

      s%extras_startup => extras_startup
      s%extras_start_step => extras_start_step
      s%extras_check_model => extras_check_model
      s%extras_finish_step => extras_finish_step
      s%extras_after_evolve => extras_after_evolve
      s%how_many_extra_history_columns => how_many_extra_history_columns
      s%data_for_extra_history_columns => data_for_extra_history_columns
      s%how_many_extra_profile_columns => how_many_extra_profile_columns
      s%data_for_extra_profile_columns => data_for_extra_profile_columns

      ! pulsation info
      s%other_photo_write => photo_write
      s%other_photo_read => photo_read

      ! this is optional
      s%other_wind => brott_wind
      s%other_adjust_mdot => my_adjust_mdot
      s%other_before_struct_burn_mix => my_before_struct_burn_mix
      s%other_kap_get => my_other_kap_get

      ! store user provided options from the inlist

      in_inlist_pulses = s%x_logical_ctrl(22)
      max_dt_before_pulse = s%x_ctrl(17)
      max_dt_during_pulse = s%x_ctrl(18)
      remesh_for_envelope_model = s%x_logical_ctrl(23)
      turn_off_remesh = s%x_logical_ctrl(24)
      kick_model_number = s%x_ctrl(11)
      timestep_drop_model_number = s%x_ctrl(13)
      turn_off_remesh_model_number = s%x_ctrl(12)
   end subroutine extras_controls

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
      type(star_info), pointer :: s
      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return

      L1 = Lsurf
      M1 = Msurf
      R1 = Rsurf
      T1 = Tsurf

      h1 = s%net_iso(ih1)
      he4 = s%net_iso(ihe4)
      Xs = s%xa(h1, 1)
      Ys = s%xa(he4, 1)
      ! Z=0.0142 is Z from Asplund et al. 2009
      Z_div_Z_solar = s%kap_rq%Zbase/0.0142d0
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
         w = alfa*vink_wind + (1d0 - alfa)*hamann_wind
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
               alfa = (T1 - (Teff_jump - dT))/(2*dT)
            end if
         end if

         if (alfa > 0) then ! eval hot side wind (eqn 24)
            vinf_div_vesc = 2.6d0 ! this is the hot side galactic value
            vinf_div_vesc = vinf_div_vesc*pow(Z_div_Z_solar, 0.13d0) ! corrected for Z
            logMdot = &
               -6.697d0 &
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
            vinf_div_vesc = vinf_div_vesc*pow(Z_div_Z_solar, 0.13d0) ! corrected for Z
            logMdot = &
               -6.688d0 &
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
                  + 1.24d0*log10(L1/Lsun) &
                  + 0.16d0*log10(M1/Msun) &
                  + 0.81d0*log10(R1/Rsun) &
                  + 0.85d0*log10(Z_div_Z_solar)
         w = exp10(log10w)
      end subroutine eval_Nieuwenhuijzen_wind

      subroutine eval_Hamann_wind(w)
         ! Hamann, W.-R.; Koesterke, L.; Wessolowski, U. 1995, A&A, 299, 151
         real(dp), intent(out) :: w
         real(dp) :: log10w
         include 'formats'
         log10w = -11.95d0 &
                  + 1.5d0*log10(L1/Lsun) &
                  - 2.85d0*Xs &
                  + 0.85d0*log10(Z_div_Z_solar)
         w = exp10(log10w)
      end subroutine eval_Hamann_wind

   end subroutine brott_wind

   subroutine my_adjust_mdot(id, ierr)
      use star_def
      integer, intent(in) :: id
      integer, intent(out) :: ierr
      type(star_info), pointer :: s
      real(dp) :: Lrad_div_Ledd
      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return
      if (s%generations > 2) then
         write (*, *) "check mdots", s%mstar_dot, s%mstar_dot_old
         if (abs(s%mstar_dot) > 1.05d0*abs(s%mstar_dot_old)) then
            s%mstar_dot = 1.05d0*s%mstar_dot_old
         else if (abs(s%mstar_dot) < 0.95d0*abs(s%mstar_dot_old)) then
            s%mstar_dot = 0.95d0*s%mstar_dot_old
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

      type(star_info), pointer :: s
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
      if (k == 1 .and. s%u_flag) then !very surface cell can go mad, things are more stable if we fix opacity
         ! this is to support restarts, as xh_start and r_start are
         ! not properly set when model loads
         if (s%solver_iter > 0) then
            velocity = s%xh_start(s%i_u, 1)
            radius = s%r_start(1)
         else
            velocity = s%xh(s%i_u, 1)
            radius = s%r(1)
         end if
         if (velocity > sqrt(2*s%cgrav(1)*s%m(1)/radius)) then
            kap = 0.2d0*(1 + s%X(1))
            dln_kap_dlnRho = 0d0
            dln_kap_dlnT = 0d0
            return
         else
            call kap_get( &
               s%kap_handle, species, chem_id, net_iso, xa, &
               log10_rho, log10_T, lnfree_e, d_lnfree_e_dlnRho, d_lnfree_e_dlnT, &
               eta, d_eta_dlnRho, d_eta_dlnT, &
               kap_fracs, kap, dln_kap_dlnRho, dln_kap_dlnT, dln_kap_dxa, ierr)
         end if
      else
         call kap_get( &
            s%kap_handle, species, chem_id, net_iso, xa, &
            log10_rho, log10_T, lnfree_e, d_lnfree_e_dlnRho, d_lnfree_e_dlnT, &
            eta, d_eta_dlnRho, d_eta_dlnT, &
            kap_fracs, kap, dln_kap_dlnRho, dln_kap_dlnT, dln_kap_dxa, ierr)
      end if

   end subroutine my_other_kap_get

   subroutine extras_startup(id, restart, ierr)
      integer, intent(in) :: id
      logical, intent(in) :: restart
      integer, intent(out) :: ierr
      type(star_info), pointer :: s
      include 'formats'
      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return
      call test_suite_startup(s, restart, ierr)
      call TDC_pulsation_extras_startup(id, restart, ierr)

      ! Initialize GYRE

      call init('gyre.in')

      ! Set constants

      call set_constant('G_GRAVITY', standard_cgrav)
      call set_constant('C_LIGHT', clight)
      call set_constant('A_RADIATION', crad)

      call set_constant('M_SUN', Msun)
      call set_constant('R_SUN', Rsun)
      call set_constant('L_SUN', Lsun)

      call set_constant('GYRE_DIR', TRIM(mesa_dir)//'/build/gyre/src')

      if (.not. restart .and. in_inlist_pulses) then
          initial_model_number = s% model_number
      end if
      !initial_model_number = 0 ! since we are setting model # to 0 in inlist_pulses

      ! for rsp style mesh
      if (.not. restart .and. in_inlist_pulses .and. remesh_for_envelope_model) then
         call remesh_for_TDC_pulsation(id, ierr)
      end if
   end subroutine extras_startup

   subroutine extras_after_evolve(id, ierr)
      integer, intent(in) :: id
      integer, intent(out) :: ierr
      type(star_info), pointer :: s
      real(dp) :: dt
      character(len=strlen) :: test
      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return
      call test_suite_after_evolve(s, ierr)

      if (.not. s%x_logical_ctrl(37)) return
      call final()
   end subroutine extras_after_evolve

   ! returns either keep_going, retry, or terminate.
   integer function extras_check_model(id)
      integer, intent(in) :: id
      integer :: ierr, k
      real(dp) :: max_v
      type(star_info), pointer :: s
      include 'formats'
      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return
      extras_check_model = keep_going

   end function extras_check_model

   integer function how_many_extra_history_columns(id)
      integer, intent(in) :: id
      integer :: ierr
      type(star_info), pointer :: s
      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return
      how_many_extra_history_columns = TDC_pulsation_how_many_extra_history_columns(id)
   end function how_many_extra_history_columns

   subroutine data_for_extra_history_columns(id, n, names, vals, ierr)
      integer, intent(in) :: id, n
      character(len=maxlen_history_column_name) :: names(n)
      real(dp) :: vals(n), v_esc
      integer, intent(out) :: ierr
      type(star_info), pointer :: s
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
      type(star_info), pointer :: s
      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return
      how_many_extra_profile_columns = TDC_pulsation_how_many_extra_profile_columns(id)

   end function how_many_extra_profile_columns

   subroutine data_for_extra_profile_columns(id, n, nz, names, vals, ierr)
      use star_def, only: star_info, maxlen_profile_column_name
      use const_def, only: dp
      integer, intent(in) :: id, n, nz
      character(len=maxlen_profile_column_name) :: names(n)
      real(dp) :: vals(nz, n)
      integer, intent(out) :: ierr
      type(star_info), pointer :: s
      integer :: k
      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return
      call TDC_pulsation_data_for_extra_profile_columns(id, n, nz, names, vals, ierr)

   end subroutine data_for_extra_profile_columns

   integer function extras_start_step(id)
      integer, intent(in) :: id
      integer :: ierr
      type(star_info), pointer :: s
      include 'formats'
      extras_start_step = terminate
      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return

      !this is used to ensure we read the right inlist options
      s%use_other_before_struct_burn_mix = .true.

      ! we want to ignore T gradient equation for a few steps after remesh
      if (s%model_number < initial_model_number + 10 .and. in_inlist_pulses) then
         s%convergence_ignore_equL_residuals = .true.
      else if (in_inlist_pulses) then
         s%convergence_ignore_equL_residuals = .false.
      end if

      if (s%model_number == kick_model_number .and. in_inlist_pulses &
          .and. s%x_logical_ctrl(5)) then

         ! if v= 0, turn on v so we can kick
         if (.not. s%v_flag .and. .not. s%u_flag) then
            call star_set_v_flag(id, .true., ierr)
         end if

         call gyre_in_mesa_extras_set_velocities(s, .false., ierr)
         write (*, *) 'kick'
         write (*, *) 'kick'
         write (*, *) 'kick'
         write (*, *) 'kick'
         write (*, *) 'kick'

      end if

      call my_before_struct_burn_mix(s%id, s%dt, extras_start_step)

      ! add stopping condition for testing.
      if ((.not. in_inlist_pulses) .and. s%center_he4 < 2d-1) then
         s%Teff_lower_limit = exp10(3.75d0)
      else
         s%Teff_lower_limit = -1d99
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
      type(star_info), pointer :: s
      include 'formats'
      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return

      if (in_inlist_pulses) then
         if (s%model_number > timestep_drop_model_number) then
            s%max_timestep = max_dt_during_pulse
         else
            s%max_timestep = max_dt_before_pulse
         end if

         ! time step control on pulsations
         if (period > 0d0 .and. period/s%max_timestep < 600 .and. &
             s%model_number > timestep_drop_model_number) then
            s%max_timestep = period/600d0
         end if

         if (s%model_number > turn_off_remesh_model_number .and. turn_off_remesh) then
            s%okay_to_remesh = .false.
         end if
      end if

      ! reading inlists can turn this flag off for some reason
      s%use_other_before_struct_burn_mix = .true.

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
      integer :: ierr, k
      real(dp) :: max_vel_inside, vesc_for_cell, vesc_surf !check_avg_v_div_vesc
      type(star_info), pointer :: s
      include 'formats'
      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return

      extras_finish_step = keep_going
      extras_finish_step = TDC_pulsation_extras_finish_step(id)

!         if (.not. s% x_logical_ctrl(37)) return
!         extras_finish_step = gyre_in_mesa_extras_finish_step(id)

      if (extras_finish_step == terminate) s%termination_code = t_extras_finish_step

   end function extras_finish_step


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

