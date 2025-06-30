! ***********************************************************************
!
!   Copyright (C) 2022  The MESA Team
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

      module run_binary
      implicit none

      contains

      subroutine do_run_binary(tst)
         use binary_lib, only: run1_binary
         use run_star_extras
         use run_binary_extras

         logical, intent(in) :: tst

         integer :: ierr

         call run1_binary(tst, &
            ! star extras
            extras_controls, &
            ! binary extras
            extras_binary_controls, &
            ierr)

      end subroutine do_run_binary

      end module run_binary

