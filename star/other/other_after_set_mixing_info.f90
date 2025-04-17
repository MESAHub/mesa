! ***********************************************************************
!
!   Copyright (C) 2014  The MESA Team
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

module other_after_set_mixing_info

   implicit none

contains

   subroutine default_other_after_set_mixing_info(id, ierr)
      use const_def, only: dp
      use star_def
      integer, intent(in) :: id
      integer, intent(out) :: ierr
      ierr = 0
   end subroutine default_other_after_set_mixing_info

end module other_after_set_mixing_info

