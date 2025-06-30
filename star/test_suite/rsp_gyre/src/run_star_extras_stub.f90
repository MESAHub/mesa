! ***********************************************************************
!
!   Copyright (C) 2018-2019  Rich Townsend & The MESA Team
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
      use gyre_lib
      use math_lib

      implicit none

      include "test_suite_extras_def.inc"

      contains

      include "test_suite_extras.inc"


      include 'run_star_extras.inc'


      subroutine extras_controls(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s

         write(*,*) 'cannot run rsp without gyre.'
         write(*,*) 'this test was intentionally skipped'
         write(*,*) 'good match for period', -1d0, -1d0

         open(unit=30, file='final.mod', action='write', status='replace')
         write(30,*) 'fake final.mod'
         close(30)

         ierr = -1
         return

      end subroutine extras_controls


      end module run_star_extras

