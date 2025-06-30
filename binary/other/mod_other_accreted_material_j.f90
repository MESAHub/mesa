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

      module mod_other_accreted_material_j

      ! NOTE: remember to set true:
      ! use_other_accreted_material_j = .true.

      implicit none

      contains

      subroutine null_other_accreted_material_j(binary_id, ierr)
         use binary_def, only : binary_info, binary_ptr
         use const_def, only: dp
         integer, intent(in) :: binary_id
         integer, intent(out) :: ierr
         type (binary_info), pointer :: b
         ierr = 0
         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in binary_ptr'
            return
         end if
         b% accretion_mode = 0
         b% acc_am_div_kep_am = 0.0d0

      end subroutine null_other_accreted_material_j

      end module mod_other_accreted_material_j

