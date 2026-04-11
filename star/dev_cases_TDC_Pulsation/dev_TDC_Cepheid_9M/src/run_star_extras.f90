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
   use gyre_mesa_m

   use interp_1d_def, only: pm_work_size
   use interp_1d_lib, only: interp_pm, interp_values, interp_value

   implicit none

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


   subroutine extras_startup(id, restart, ierr)
      integer, intent(in) :: id
      logical, intent(in) :: restart
      integer, intent(out) :: ierr
      type(star_info), pointer :: s
      include 'formats'
      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return
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
      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return

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
            write (*,*) 'turning on v_flag hydro for kick'
         end if

         call gyre_in_mesa_extras_set_velocities(s, .false., ierr)
         write (*, *) 'kick'
         write (*, *) 'kick'
         write (*, *) 'kick'
         write (*, *) 'kick'
         write (*, *) 'kick'

      end if

      call my_before_struct_burn_mix(s%id, s%dt, extras_start_step)

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
