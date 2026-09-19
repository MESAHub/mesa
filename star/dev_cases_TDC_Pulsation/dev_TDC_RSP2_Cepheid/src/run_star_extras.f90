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

   integer, parameter :: num_rsp2_profile_columns = 14

   include 'run_star_extras_TDC_pulsation_defs.inc'

   logical :: in_inlist_pulses, turn_off_remesh
   integer :: steps_per_period, timestep_drop_model_number, &
      turn_off_remesh_model_number
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

      in_inlist_pulses = s% x_logical_ctrl(22)
      turn_off_remesh = s% x_logical_ctrl(24)
      steps_per_period = s% x_integer_ctrl(8)
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
      integer :: ierr
      type(star_info), pointer :: s

      how_many_extra_profile_columns = &
         TDC_pulsation_how_many_extra_profile_columns(id)
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return
      if (s% RSP2_flag) how_many_extra_profile_columns = &
         how_many_extra_profile_columns + num_rsp2_profile_columns
   end function how_many_extra_profile_columns


   subroutine data_for_extra_profile_columns(id, n, nz, names, vals, ierr)
      integer, intent(in) :: id, n, nz
      character(len=maxlen_profile_column_name) :: names(n)
      real(dp) :: vals(nz, n)
      integer, intent(out) :: ierr
      type(star_info), pointer :: s
      integer :: i, j, k, km1, num_tdc_columns
      real(dp) :: alfa, beta, PII_div_Hp

      num_tdc_columns = TDC_pulsation_how_many_extra_profile_columns(id)
      call TDC_pulsation_data_for_extra_profile_columns( &
         id, num_tdc_columns, nz, names(1:num_tdc_columns), vals(:,1:num_tdc_columns), ierr)
      if (ierr /= 0) return
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return
      if (.not. s% RSP2_flag) return

      i = num_tdc_columns + 1
      names(i) = 'w_start'; vals(:,i) = s% w_start(1:nz); i = i+1
      names(i) = 'Y_face_start'; vals(:,i) = s% Y_face_start(1:nz); i = i+1
      names(i) = 'L_start'; vals(:,i) = s% L_start(1:nz); i = i+1
      names(i) = 'Lt_start'; vals(:,i) = s% Lt_start(1:nz); i = i+1
      names(i) = 'rho_start'; vals(:,i) = s% rho_start(1:nz); i = i+1
      names(i) = 'energy_start'; vals(:,i) = s% energy_start(1:nz); i = i+1
      names(i) = 'r_start'; vals(:,i) = s% r_start(1:nz); i = i+1
      names(i) = 'v_start'; vals(:,i) = 0d0
      if (s% v_flag) vals(:,i) = s% v_start(1:nz)
      i = i+1
      names(i) = 'u_start'; vals(:,i) = 0d0
      if (s% u_flag) vals(:,i) = s% u_start(1:nz)
      i = i+1
      names(i) = 'T_start'; vals(:,i) = s% T_start(1:nz); i = i+1
      names(i) = 'csound_start'; vals(:,i) = s% csound_start(1:nz); i = i+1

      j = i
      names(j) = 'w_face'
      names(j+1) = 'Source_div_w'
      names(j+2) = 'Ptrb_rsp2'
      vals(:,j:j+2) = 0d0
      do k = 1, nz
         km1 = max(1,k-1)
         if (k == 1) then
            alfa = 1d0
         else if (s% RSP2_use_mass_interp_face_values) then
            alfa = s% dq(k-1)/(s% dq(k-1) + s% dq(k))
         else
            alfa = 0.5d0
         end if
         beta = 1d0 - alfa
         vals(k,j) = alfa*s% w(k) + beta*s% w(km1)
         if (s% mixing_length_alpha == 0d0 .or. &
               k <= s% RSP2_num_outermost_cells_forced_nonturbulent .or. &
               k > nz - int(nz/s% RSP2_nz_div_IBOTOM)) cycle

         ! Match compute_Source_div_w without dividing by a possibly zero w.
         PII_div_Hp = s% PII(k)/s% Hp_face(k)
         if (k < nz) PII_div_Hp = 0.5d0*(PII_div_Hp + s% PII(k+1)/s% Hp_face(k+1))
         vals(k,j+1) = PII_div_Hp*s% Peos(k)*s% chiT(k)/(s% rho(k)*s% chiRho(k)*s% Cp(k))
         vals(k,j+2) = s% RSP2_alfap*(2d0/3d0)*s% rho(k)*s% w(k)**2
      end do
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

end module run_star_extras
