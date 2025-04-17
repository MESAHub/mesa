! ***********************************************************************
!
!   Copyright (C) 2012  The MESA Team
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

module other_adjust_mdot

   ! consult star/other/README for general usage instructions
   ! control name: use_other_adjust_mdot = .true.
   ! procedure pointer: s% other_adjust_mdot => my_routine

   implicit none

contains

   ! your routine will be called after winds and before mass adjustment
   subroutine null_other_adjust_mdot(id, ierr)
      use star_def
      integer, intent(in) :: id
      integer, intent(out) :: ierr
      ierr = 0
   end subroutine null_other_adjust_mdot

end module other_adjust_mdot

