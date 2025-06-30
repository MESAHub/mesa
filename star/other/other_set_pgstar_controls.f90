! ***********************************************************************
!
!   Copyright (C) 2021  The MESA Team
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

module other_set_pgstar_controls

   ! consult star/other/README for general usage instructions
   ! control name: use_other_set_pgstar_controls = .true.
   ! procedure pointer: s% other_set_pgstar_controls => my_routine

   implicit none

contains

   subroutine default_other_set_pgstar_controls(id)
      use star_def
      integer, intent(in) :: id
      type(star_info), pointer :: s
      integer :: ierr
      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return
   end subroutine default_other_set_pgstar_controls

end module other_set_pgstar_controls

