! ***********************************************************************
!
!   Copyright (C) 2012  Bill Paxton, Pablo Marchant & The MESA Team
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

      module mod_other_binary_ce

      ! NOTE: remember to set true:
      ! use_other_CE_rlo_mdot = .true.

      implicit none

      contains

      subroutine null_other_CE_init(binary_id, restart, ierr)
         use const_def, only: dp
         use star_def
         use binary_def, only : binary_info, binary_ptr
         integer, intent(in) :: binary_id
         logical, intent(in) :: restart
         integer, intent(out) :: ierr
         type (binary_info), pointer :: b
         write(*,*) "WARNING: using null_other_CE_init"
         ierr = 0
         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in binary_ptr'
            return
         end if
      end subroutine null_other_CE_init

      subroutine null_other_CE_rlo_mdot(binary_id, mdot, ierr)
         use const_def, only: dp
         use star_def
         use binary_def, only : binary_info, binary_ptr
         integer, intent(in) :: binary_id
         real(dp), intent(out) :: mdot
         integer, intent(out) :: ierr
         type (binary_info), pointer :: b
         write(*,*) "WARNING: using null_other_CE_rlo_mdot"
         ierr = 0
         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in binary_ptr'
            return
         end if
         mdot = -1d-99
      end subroutine null_other_CE_rlo_mdot

      integer function null_other_CE_binary_evolve_step(binary_id)
         use const_def, only: dp
         use star_def
         use binary_def, only : binary_info, binary_ptr
         integer, intent(in) :: binary_id
         integer :: ierr
         type (binary_info), pointer :: b
         write(*,*) "WARNING: using null_other_CE_binary_evolve_step"
         ierr = 0
         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in binary_ptr'
            return
         end if
         null_other_CE_binary_evolve_step = keep_going
      end function null_other_CE_binary_evolve_step

      integer function null_other_CE_binary_finish_step(binary_id)
         use star_def
         use binary_def, only : binary_info, binary_ptr
         integer, intent(in) :: binary_id
         integer :: ierr
         type (binary_info), pointer :: b
         write(*,*) "WARNING: using null_other_CE_binary_check_model"
         ierr = 0
         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in binary_ptr'
            return
         end if
         null_other_CE_binary_finish_step = keep_going
      end function null_other_CE_binary_finish_step

      end module mod_other_binary_ce

