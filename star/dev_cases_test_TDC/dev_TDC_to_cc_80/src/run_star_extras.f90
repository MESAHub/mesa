! ***********************************************************************
!
!   Copyright (C) 2010  Bill Paxton & The MESA Team
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
      use gyre_lib

      implicit none

      include "test_suite_extras_def.inc"

! here are the x controls used below

!alpha_mlt_routine
         !alpha_H = s% x_ctrl(21)
         !alpha_other = s% x_ctrl(22)
         !H_limit = s% x_ctrl(23)

!gyre
      !x_logical_ctrl(37) = .false. ! if true, then run GYRE
      !x_integer_ctrl(1) = 2 ! output GYRE info at this step interval
      !x_logical_ctrl(1) = .false. ! save GYRE info whenever save profile
      !x_integer_ctrl(2) = 2 ! max number of modes to output per call
      !x_logical_ctrl(2) = .false. ! output eigenfunction files
      !x_integer_ctrl(3) = 0 ! mode l (e.g. 0 for p modes, 1 for g modes)
      !x_integer_ctrl(4) = 1 ! order
      !x_ctrl(1) = 0.158d-05 ! freq ~ this (Hz)
      !x_ctrl(2) = 0.33d+03 ! growth < this (days)


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
         s% other_alpha_mlt => alpha_mlt_routine
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


      subroutine extras_startup(id, restart, ierr)
         integer, intent(in) :: id
         logical, intent(in) :: restart
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         call test_suite_startup(s, restart, ierr)

         if (.not. s% x_logical_ctrl(37)) return

         ! Initialize GYRE

         call gyre_init('gyre.in')

         ! Set constants

         call gyre_set_constant('G_GRAVITY', standard_cgrav)
         call gyre_set_constant('C_LIGHT', clight)
         call gyre_set_constant('A_RADIATION', crad)

         call gyre_set_constant('M_SUN', Msun)
         call gyre_set_constant('R_SUN', Rsun)
         call gyre_set_constant('L_SUN', Lsun)

         call gyre_set_constant('GYRE_DIR', TRIM(mesa_dir)//'/gyre/gyre')

      end subroutine extras_startup


      subroutine extras_after_evolve(id, ierr)
         use num_lib, only: find0
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         real(dp) :: dt, m
         integer :: k, nz
         include 'formats'
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         nz = s% nz
         write(*,'(A)')
         call test_suite_after_evolve(s, ierr)
         if (.not. s% x_logical_ctrl(37)) return
         call gyre_final()
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
         how_many_extra_profile_columns = 6
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

         names(1) = 'std_lg_absY'
         names(2) = 'std_gradT'
         names(3) = 'std_vc'
         names(4) = 'S0'
         names(5) = 'D0'
         names(6) = 'DR0'

         do k=1,nz
            vals(k,1) = s% xtra1_array(k)
            vals(k,2) = s% xtra2_array(k)
            vals(k,3) = s% xtra3_array(k)
            vals(k,4) = s% xtra4_array(k)
            vals(k,5) = s% xtra5_array(k)
            vals(k,6) = s% xtra6_array(k)
         end do

      end subroutine data_for_extra_profile_columns

      include 'gyre_in_mesa_extras_finish_step.inc'

      ! returns either keep_going or terminate.
      integer function extras_finish_step(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_finish_step = keep_going
         if (.not. s% x_logical_ctrl(37)) return
         extras_finish_step = gyre_in_mesa_extras_finish_step(id)
         if (extras_finish_step == terminate) &
             s% termination_code = t_extras_finish_step
      end function extras_finish_step

      end module run_star_extras

