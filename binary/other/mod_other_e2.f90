! ***********************************************************************
!
!   Copyright (C) 2022 The MESA Team
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

      module mod_other_e2

      ! NOTE: remember to set true:
      ! use_other_e2 = .true.

      implicit none


      contains

      subroutine null_other_e2(id, e2, ierr)
         use binary_def, only : binary_info, binary_ptr
         use star_def, only : star_info, star_ptr
         use const_def, only: dp
         integer, intent(in) :: id
         real(dp), intent(out) :: e2
         integer, intent(out) :: ierr
         type (binary_info), pointer :: b
         type (star_info), pointer :: s

         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in star_ptr'
            return
         end if

         call binary_ptr(s% binary_id, b, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in binary_ptr'
            return
         end if
         e2 = -1
      end subroutine null_other_e2

      end module mod_other_e2

