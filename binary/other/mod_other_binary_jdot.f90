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

      module mod_other_binary_jdot

      ! NOTE: remember to set one of:
      ! use_other_jdot_mb = .true.
      ! use_other_jdot_gr = .true.
      ! use_other_jdot_ml = .true.
      ! use_other_extra_jdot = .true.

      implicit none

      contains

      subroutine null_other_jdot_mb(binary_id, ierr)
         use binary_def, only : binary_info, binary_ptr
         integer, intent(in) :: binary_id
         integer, intent(out) :: ierr
         type (binary_info), pointer :: b
         ierr = 0
         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in binary_ptr'
            return
         end if
         b% jdot_mb = 0
      end subroutine null_other_jdot_mb

      subroutine null_other_jdot_gr(binary_id, ierr)
         use binary_def, only : binary_info, binary_ptr
         integer, intent(in) :: binary_id
         integer, intent(out) :: ierr
         type (binary_info), pointer :: b
         ierr = 0
         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in binary_ptr'
            return
         end if
         b% jdot_gr = 0
      end subroutine null_other_jdot_gr

      subroutine null_other_jdot_ml(binary_id, ierr)
         use binary_def, only : binary_info, binary_ptr
         integer, intent(in) :: binary_id
         integer, intent(out) :: ierr
         type (binary_info), pointer :: b
         ierr = 0
         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in binary_ptr'
            return
         end if
         b% jdot_ml = 0
      end subroutine null_other_jdot_ml

      subroutine null_other_extra_jdot(binary_id, ierr)
         use binary_def, only : binary_info, binary_ptr
         integer, intent(in) :: binary_id
         integer, intent(out) :: ierr
         type (binary_info), pointer :: b
         ierr = 0
         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in binary_ptr'
            return
         end if
         b% extra_jdot = 0
      end subroutine null_other_extra_jdot

      subroutine null_other_jdot_ls(binary_id, ierr)
         use binary_def, only : binary_info, binary_ptr
         integer, intent(in) :: binary_id
         integer, intent(out) :: ierr
         type (binary_info), pointer :: b
         ierr = 0
         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in binary_ptr'
            return
         end if
         b% jdot_ls = 0
      end subroutine null_other_jdot_ls

      subroutine null_other_jdot_missing_wind(binary_id, ierr)
         use binary_def, only : binary_info, binary_ptr
         integer, intent(in) :: binary_id
         integer, intent(out) :: ierr
         type (binary_info), pointer :: b
         ierr = 0
         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in binary_ptr'
            return
         end if
         b% jdot_missing_wind = 0
      end subroutine null_other_jdot_missing_wind


      end module mod_other_binary_jdot

