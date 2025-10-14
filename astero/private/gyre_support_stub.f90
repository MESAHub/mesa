! ***********************************************************************
!
!   Copyright (C) 2013  The MESA Team
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

      module gyre_support

      use star_lib
      use star_def

      implicit none

      logical, parameter :: GYRE_IS_ENABLED = .false.

      contains

      subroutine init_gyre(gyre_file,ierr)
         character(*), intent(in) :: gyre_file
         integer, intent(out) :: ierr
         ierr = -1
      end subroutine init_gyre


      subroutine do_gyre_get_modes (s, el, store_model, ierr)

         type (star_info), pointer :: s
         integer, intent(in)       :: el
         logical, intent(in)       :: store_model
         integer, intent(out)      :: ierr

         ierr = -1
      end subroutine do_gyre_get_modes


      end module gyre_support
