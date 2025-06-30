! ***********************************************************************
!
!   Copyright (C) 2010  The MESA Team
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

module other_remove_surface

   ! consult star/other/README for general usage instructions
   ! control name: use_other_remove_surface = .true.
   ! procedure pointer: s% other_remove_surface => my_routine

   implicit none

contains

   subroutine default_other_remove_surface(id, ierr, k)
      use star_def
      integer, intent(in) :: id
      integer, intent(out) :: ierr
      integer, intent(out) :: k

      type(star_info), pointer :: s
      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return

      k = 0  ! The cell to remove down to.
   end subroutine default_other_remove_surface

end module other_remove_surface

