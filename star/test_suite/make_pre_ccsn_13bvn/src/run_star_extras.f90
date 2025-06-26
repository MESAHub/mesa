! ***********************************************************************
!
!   Copyright (C) 2023  Ebraheem Farag & The MESA Team
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

      implicit none

      include "test_suite_extras_def.inc"

      integer, parameter :: I_INLIST_PART = 1  ! inlist part number

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


      !integer function extras_startup(id, restart, ierr)
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


      ! returns either keep_going, retry, backup, or terminate.
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
         how_many_extra_history_columns = 1
      end function how_many_extra_history_columns


      subroutine data_for_extra_history_columns(id, n, names, vals, ierr)
         integer, intent(in) :: id, n
         character (len=maxlen_history_column_name) :: names(n)
         real(dp) :: vals(n)
         integer :: k
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         names(1) = 'log_total_drag_energy'
         vals(1) = 0
         if (s% v_flag .and. s% drag_coefficient > 0) then
            do k=1,s% nz
               if (s% q(k) >= s% min_q_for_drag) then
                  vals(1) = s% FdotV_drag_energy(k) * s% dm(k) + vals(1)
               else
                  vals(1) = vals(1)
               end if
            end do
            vals(1) = log10(max(1d-40,abs(vals(1))))
         else
            vals(1) = 0d0
         end if

      end subroutine data_for_extra_history_columns


      integer function how_many_extra_profile_columns(id)
         use star_def, only: star_info
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_profile_columns = 3
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
         names(1) = 'zbar_div_abar'
         names(2) = 'dvdt_drag'
         names(3) = 'FdotV_drag_energy'

         if(s% v_flag .and. s% drag_coefficient > 0) then
            do k=1,s% nz
               vals(k,1) = s% zbar(k)/s% abar(k)
               vals(k,2) = s% dvdt_drag(k)
               vals(k,3) = s% FdotV_drag_energy(k)
            end do
         else
            do k=1,s% nz
               vals(k,1) = s% zbar(k)/s% abar(k)
               vals(k,2) = 0d0
               vals(k,3) = 0d0
            end do
         end if

      end subroutine data_for_extra_profile_columns

 ! returns either keep_going or terminate.
      integer function extras_start_step(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s

         include 'formats'
         extras_start_step = keep_going
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

      end function extras_start_step


      ! returns either keep_going or terminate.
      integer function extras_finish_step(id)
         integer, intent(in) :: id
         integer :: ierr
         real(dp) :: env_mass_check
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_finish_step = keep_going
         ! custom mass_change for inlist_remove written by -EbF
         env_mass_check = s% star_mass - s% he_core_mass
         if (s% x_logical_ctrl(1)) then
             if (env_mass_check >= 1d-1 ) then
                 s% mass_change = -1d-3
             else if (env_mass_check <= 1d-1 .and. env_mass_check >= s% star_species_mass_min_limit) then
                s% mass_change = -1d-4
             end if
         else
             s% mass_change = 0d0
         end if
      end function extras_finish_step

      subroutine extras_photo_read(id, iounit, ierr)
         integer, intent(in) :: id, iounit
         integer, intent(out) :: ierr
         integer :: inlist_part
         type (star_info), pointer :: s
         ierr = 0

         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         read(iounit,iostat=ierr)  inlist_part

         if(inlist_part/= s% x_integer_ctrl(I_INLIST_PART)) then
            call mesa_error(__FILE__,__LINE__,'Error: Photo was saved for different inlist')
            ierr=-1
            return
         end if

      end subroutine extras_photo_read

      subroutine extras_photo_write(id, iounit)
         integer, intent(in) :: id, iounit
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0

         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         write(iounit) s% x_integer_ctrl(I_INLIST_PART)

      end subroutine extras_photo_write

   end module run_star_extras

