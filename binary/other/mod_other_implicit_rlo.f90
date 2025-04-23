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

      module mod_other_implicit_rlo

      ! NOTE: remember to set one of these to true:
      ! use_other_check_implicit_rlo_mdot = .true.
      ! use_other_implicit_function_to_solve = .true.

      use const_def, only: dp

      implicit none

      contains

      integer function null_other_check_implicit_rlo(binary_id, new_mdot)
         use binary_def, only : binary_info, binary_ptr
         use const_def, only: dp
         use star_def
         integer, intent(in) :: binary_id
         real(dp), intent(out) :: new_mdot
         integer :: ierr
         type (binary_info), pointer :: b
         ierr = 0
         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in binary_ptr'
            return
         end if
         new_mdot = 0d0
         null_other_check_implicit_rlo = keep_going
         write(*,*) "WARNING: using null_other_check_implicit_rlo"
      end function null_other_check_implicit_rlo

      subroutine null_other_implicit_function_to_solve(binary_id, &
         function_to_solve, use_sum, detached, ierr)
         use binary_def, only : binary_info, binary_ptr
         use const_def, only: dp
         integer, intent(in) :: binary_id
         real(dp), intent(out) :: function_to_solve
         logical, intent(out) :: use_sum, detached
         integer, intent(out) :: ierr
         type (binary_info), pointer :: b
         ierr = 0
         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in binary_ptr'
            return
         end if
         function_to_solve = 0d0
         use_sum = .false.
         write(*,*) "WARNING: using null_other_implicit_function_to_solve"
      end subroutine null_other_implicit_function_to_solve

      end module mod_other_implicit_rlo

