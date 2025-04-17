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

module other_brunt

   ! consult star/other/README for general usage instructions
   ! control name: use_other_brunt = .true.
   ! procedure pointer: s% other_brunt => my_routine

   implicit none

contains

   subroutine default_other_brunt(id, ierr)
      use star_def
      integer, intent(in) :: id
      integer, intent(out) :: ierr
      type(star_info), pointer :: s
      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return
      s%brunt_N2(1:s%nz) = 0
   end subroutine default_other_brunt

end module other_brunt

