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

module other_cgrav

   ! consult star/other/README for general usage instructions
   ! control name: use_other_cgrav = .true.
   ! procedure pointer: s% other_cgrav => my_routine

   ! note that other_cgrav only changes G in the stellar structure
   ! the binary module is unaffected by changes in cgrav

   implicit none

contains

   subroutine default_other_cgrav(id, ierr)
      use star_def
      use const_def, only: standard_cgrav
      integer, intent(in) :: id
      integer, intent(out) :: ierr
      type(star_info), pointer :: s
      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return
      s%cgrav(:) = standard_cgrav
   end subroutine default_other_cgrav

end module other_cgrav
