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
   use utils_lib, only: is_bad

   implicit none

   include 'run_star_extras_TDC_pulsation_defs.inc'

   logical :: in_inlist_pulses, turn_off_remesh
   integer :: steps_per_period, timestep_drop_model_number, &
      turn_off_remesh_model_number, freeze_surface_cell_remesh_model_number
   real(dp) :: max_dt_before_pulse, max_dt_during_pulse

contains

   include 'run_star_extras_TDC_pulsation.inc'

   subroutine extras_controls(id, ierr)
      integer, intent(in) :: id
      integer, intent(out) :: ierr
      type(star_info), pointer :: s

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
      s% other_photo_write => photo_write
      s% other_photo_read => photo_read
      s% other_kap_get => my_other_kap_get

      in_inlist_pulses = s% x_logical_ctrl(22)
      turn_off_remesh = s% x_logical_ctrl(24)
      steps_per_period = s% x_integer_ctrl(8)
      freeze_surface_cell_remesh_model_number = s% x_integer_ctrl(9)
      timestep_drop_model_number = int(s% x_ctrl(13))
      turn_off_remesh_model_number = int(s% x_ctrl(12))
      max_dt_before_pulse = s% x_ctrl(17)
      max_dt_during_pulse = s% x_ctrl(18)
   end subroutine extras_controls


   subroutine extras_startup(id, restart, ierr)
      integer, intent(in) :: id
      logical, intent(in) :: restart
      integer, intent(out) :: ierr

      call TDC_pulsation_extras_startup(id, restart, ierr)
   end subroutine extras_startup


   integer function extras_start_step(id)
      integer, intent(in) :: id
      integer :: ierr
      real(dp) :: dt_limit
      logical :: have_dt_limit
      type(star_info), pointer :: s

      extras_start_step = terminate
      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return

      if (in_inlist_pulses .and. s% x_logical_ctrl(25)) &
         s% convergence_ignore_equL_residuals = s% model_number < 10

      if (in_inlist_pulses) then
         if (s% model_number > timestep_drop_model_number) then
            if (max_dt_during_pulse > 0d0) &
               s% max_timestep = max_dt_during_pulse
         else
            if (max_dt_before_pulse > 0d0) &
               s% max_timestep = max_dt_before_pulse
         end if

         have_dt_limit = .false.
         dt_limit = 0d0
         if (steps_per_period > 0 .and. &
               s% model_number > timestep_drop_model_number) then
            if (num_periods < 1) then
               if (.not. is_bad(s% dynamic_timescale) .and. &
                     s% dynamic_timescale > 0d0) then
                  dt_limit = s% dynamic_timescale/real(steps_per_period, dp)
                  have_dt_limit = .true.
               end if
            else if (period > 0d0) then
               dt_limit = period/real(steps_per_period, dp)
               have_dt_limit = .true.
            end if
         end if

         if (have_dt_limit) then
            if (s% max_timestep <= 0d0) then
               s% max_timestep = dt_limit
            else
               s% max_timestep = min(s% max_timestep, dt_limit)
            end if
         end if

         if (s% model_number > turn_off_remesh_model_number .and. &
               turn_off_remesh) s% okay_to_remesh = .false.

         if (freeze_surface_cell_remesh_model_number >= 0 .and. &
               s% model_number >= freeze_surface_cell_remesh_model_number) then
            s% split_merge_amr_okay_to_split_1 = .false.
            s% merge_amr_ignore_surface_cells = .true.
         end if
      end if

      extras_start_step = keep_going
   end function extras_start_step


   integer function extras_check_model(id)
      integer, intent(in) :: id
      integer :: ierr
      type(star_info), pointer :: s

      extras_check_model = terminate
      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return

      extras_check_model = keep_going
   end function extras_check_model


   integer function extras_finish_step(id)
      integer, intent(in) :: id
      integer :: ierr
      type(star_info), pointer :: s

      extras_finish_step = terminate
      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return

      extras_finish_step = TDC_pulsation_extras_finish_step(id)
      if (extras_finish_step == terminate) &
         s% termination_code = t_extras_finish_step
   end function extras_finish_step


   subroutine extras_after_evolve(id, ierr)
      integer, intent(in) :: id
      integer, intent(out) :: ierr

      call TDC_pulsation_extras_after_evolve(id, ierr)
   end subroutine extras_after_evolve


   integer function how_many_extra_history_columns(id)
      integer, intent(in) :: id

      how_many_extra_history_columns = &
         TDC_pulsation_how_many_extra_history_columns(id)
   end function how_many_extra_history_columns


   subroutine data_for_extra_history_columns(id, n, names, vals, ierr)
      integer, intent(in) :: id, n
      character(len=maxlen_history_column_name) :: names(n)
      real(dp) :: vals(n)
      integer, intent(out) :: ierr

      call TDC_pulsation_data_for_extra_history_columns( &
         id, n, names, vals, ierr)
   end subroutine data_for_extra_history_columns


   integer function how_many_extra_profile_columns(id)
      integer, intent(in) :: id

      how_many_extra_profile_columns = &
         TDC_pulsation_how_many_extra_profile_columns(id)
   end function how_many_extra_profile_columns


   subroutine data_for_extra_profile_columns(id, n, nz, names, vals, ierr)
      integer, intent(in) :: id, n, nz
      character(len=maxlen_profile_column_name) :: names(n)
      real(dp) :: vals(nz, n)
      integer, intent(out) :: ierr

      call TDC_pulsation_data_for_extra_profile_columns( &
         id, n, nz, names, vals, ierr)
   end subroutine data_for_extra_profile_columns


   subroutine photo_write(id, iounit)
      integer, intent(in) :: id, iounit

      call TDC_pulsation_photo_write(id, iounit)
   end subroutine photo_write


   subroutine photo_read(id, iounit, ierr)
      integer, intent(in) :: id, iounit
      integer, intent(out) :: ierr

      call TDC_pulsation_photo_read(id, iounit, ierr)
   end subroutine photo_read


   subroutine my_other_kap_get( &
         id, k, handle, species, chem_id, net_iso, xa, &
         log10_rho, log10_T, &
         lnfree_e, d_lnfree_e_dlnRho, d_lnfree_e_dlnT, &
         eta, d_eta_dlnRho, d_eta_dlnT, &
         kap_fracs, kap, dln_kap_dlnRho, dln_kap_dlnT, dln_kap_dxa, ierr)

      use kap_def, only: num_kap_fracs
      use kap_lib
      use chem_lib, only: basic_composition_info
      use eos_def, only: i_frac_ideal

      real(dp), parameter :: kap_logT_floor = 2.7001d0
      real(dp), parameter :: kap_logT_blend_width = 0.1d0
      real(dp), parameter :: kap_logT_blend_top = &
         kap_logT_floor + kap_logT_blend_width

      ! INPUT
      integer, intent(in) :: id  ! star id if available; 0 otherwise
      integer, intent(in) :: k  ! cell number or 0 if not for a particular cell
      integer, intent(in) :: handle  ! from alloc_kap_handle
      integer, intent(in) :: species
      integer, pointer :: chem_id(:)  ! maps species to chem id
         ! index from 1 to species
         ! value is between 1 and num_chem_isos
      integer, pointer :: net_iso(:)  ! maps chem id to species number
         ! index from 1 to num_chem_isos (defined in chem_def)
         ! value is 0 if the iso is not in the current net
         ! else is value between 1 and number of species in current net
      real(dp), intent(in) :: xa(:)  ! mass fractions
      real(dp), intent(in) :: log10_rho  ! density
      real(dp), intent(in) :: log10_T  ! temperature
      real(dp), intent(in) :: lnfree_e, d_lnfree_e_dlnRho, d_lnfree_e_dlnT
         ! free_e := total combined number per nucleon of free electrons and positrons
      real(dp), intent(in) :: eta, d_eta_dlnRho, d_eta_dlnT
         ! eta := electron degeneracy parameter

      ! OUTPUT
      real(dp), intent(out) :: kap_fracs(num_kap_fracs)
      real(dp), intent(out) :: kap  ! opacity
      real(dp), intent(out) :: dln_kap_dlnRho  ! partial derivative at constant T
      real(dp), intent(out) :: dln_kap_dlnT   ! partial derivative at constant Rho
      real(dp), intent(out) :: dln_kap_dxa(:)  ! partial derivative w.r.t. to species
      integer, intent(out) :: ierr  ! 0 means AOK.

      type (star_info), pointer :: s
      real(dp) :: X, Y, Z, abar, zbar, z2bar, z53bar, ye, mass_correction, sumx
      real(dp) :: frac_ideal, lnfree_e_ideal, delta_lnfree_e_ideal
      real(dp) :: lnfree_e_for_kap
      real(dp) :: d_lnfree_e_for_kap_dlnRho, d_lnfree_e_for_kap_dlnT
      real(dp) :: kap_at_floor, dln_kap_at_floor_dlnRho
      real(dp) :: ln_kap, ln_kap_at_floor
      real(dp) :: blend, dblend_dlnT, blend_coordinate
      real(dp) :: kap_fracs_at_floor(num_kap_fracs)
      real(dp) :: dln_kap_dxa_at_floor(size(dln_kap_dxa))

      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return

      kap_fracs = 0
      kap = 0; dln_kap_dlnRho = 0; dln_kap_dlnT = 0; dln_kap_dxa = 0

      lnfree_e_for_kap = lnfree_e
      d_lnfree_e_for_kap_dlnRho = d_lnfree_e_dlnRho
      d_lnfree_e_for_kap_dlnT = d_lnfree_e_dlnT

      if (k >= 1 .and. k <= s% nz .and. &
            log10_T > 3.7d0 .and. s% eos_frac_ideal(k) > 0d0) then
         call basic_composition_info( &
            species, chem_id, xa, X, Y, Z, abar, zbar, z2bar, z53bar, &
            ye, mass_correction, sumx)
         if (ye > 0d0) then
            frac_ideal = s% eos_frac_ideal(k)
            lnfree_e_ideal = log(1d-20/avo) - ln10*log10_rho
            delta_lnfree_e_ideal = log(ye) - lnfree_e_ideal

            ! Replace only the ideal EOS contribution to lnfree_e.
            lnfree_e_for_kap = lnfree_e + frac_ideal*delta_lnfree_e_ideal
            d_lnfree_e_for_kap_dlnRho = d_lnfree_e_dlnRho + &
               s% d_eos_dlnd(i_frac_ideal,k)*delta_lnfree_e_ideal + frac_ideal
            d_lnfree_e_for_kap_dlnT = d_lnfree_e_dlnT + &
               s% d_eos_dlnT(i_frac_ideal,k)*delta_lnfree_e_ideal
         end if
      end if

      ! The PPISN escaping-surface-cell electron-scattering override is
      ! intentionally disabled for this case.
      !if (k == 1 .and. s% u_flag) then
      !   if (s% solver_iter > 0) then
      !      velocity = s% xh_start(s% i_u,1)
      !      radius = s% r_start(1)
      !   else
      !      velocity = s% xh(s% i_u,1)
      !      radius = s% r(1)
      !   end if
      !   if (velocity > sqrt(2*s% cgrav(1)*s% m(1)/radius)) then
      !      kap = 0.2d0*(1 + s% X(1))
      !      dln_kap_dlnRho = 0d0
      !      dln_kap_dlnT = 0d0
      !      return
      !   end if
      !end if

      if (log10_T < kap_logT_blend_top) then
         ! Continue the opacity from the lowest valid table temperature.
         call kap_get( &
            handle, species, chem_id, net_iso, xa, &
            log10_rho, kap_logT_floor, lnfree_e_for_kap, &
            d_lnfree_e_for_kap_dlnRho, d_lnfree_e_for_kap_dlnT, &
            eta, d_eta_dlnRho, d_eta_dlnT, &
            kap_fracs, kap, dln_kap_dlnRho, dln_kap_dlnT, dln_kap_dxa, ierr)
         if (ierr /= 0) return

         if (log10_T <= kap_logT_floor) then
            dln_kap_dlnT = 0d0
            return
         end if

         kap_at_floor = kap
         dln_kap_at_floor_dlnRho = dln_kap_dlnRho
         kap_fracs_at_floor = kap_fracs
         dln_kap_dxa_at_floor = dln_kap_dxa

         call kap_get( &
            handle, species, chem_id, net_iso, xa, &
            log10_rho, log10_T, lnfree_e_for_kap, &
            d_lnfree_e_for_kap_dlnRho, d_lnfree_e_for_kap_dlnT, &
            eta, d_eta_dlnRho, d_eta_dlnT, &
            kap_fracs, kap, dln_kap_dlnRho, dln_kap_dlnT, dln_kap_dxa, ierr)
         if (ierr /= 0) return

         blend_coordinate = &
            (log10_T - kap_logT_floor)/kap_logT_blend_width
         blend = pow2(blend_coordinate)*(3d0 - 2d0*blend_coordinate)
         dblend_dlnT = 6d0*blend_coordinate*(1d0 - blend_coordinate)/ &
            (ln10*kap_logT_blend_width)

         ln_kap_at_floor = log(kap_at_floor)
         ln_kap = log(kap)
         dln_kap_dlnRho = (1d0 - blend)*dln_kap_at_floor_dlnRho + &
            blend*dln_kap_dlnRho
         dln_kap_dlnT = blend*dln_kap_dlnT + &
            dblend_dlnT*(ln_kap - ln_kap_at_floor)
         dln_kap_dxa = (1d0 - blend)*dln_kap_dxa_at_floor + &
            blend*dln_kap_dxa
         kap_fracs = (1d0 - blend)*kap_fracs_at_floor + blend*kap_fracs
         kap = exp((1d0 - blend)*ln_kap_at_floor + blend*ln_kap)
         return
      end if

      call kap_get( &
         handle, species, chem_id, net_iso, xa, &
         log10_rho, log10_T, lnfree_e_for_kap, &
         d_lnfree_e_for_kap_dlnRho, d_lnfree_e_for_kap_dlnT, &
         eta, d_eta_dlnRho, d_eta_dlnT, &
         kap_fracs, kap, dln_kap_dlnRho, dln_kap_dlnT, dln_kap_dxa, ierr)

   end subroutine my_other_kap_get

end module run_star_extras
